#include <cassert>
#include <memory>

#include "exp.h"
#include "scanner.h"
#include "token.h"
#include "yaml-cpp/exceptions.h"  // IWYU pragma: keep

namespace YAML {
Scanner::Scanner(std::istream& in)
    : INPUT(in),
      m_startedStream(false),
      m_endedStream(false),
      m_simpleKeyAllowed(false),
      m_canBeJSONFlow(false) {}

Scanner::~Scanner() {}

// empty
// . Returns true if there are no more tokens to be read
bool Scanner::empty() {
  EnsureTokensInQueue();
  return m_tokens.empty();
}

// pop
// . Simply removes the next token on the queue.
void Scanner::pop() {
  EnsureTokensInQueue();
  if (!m_tokens.empty())
    m_tokens.pop();
}

// peek
// . Returns (but does not remove) the next token on the queue.
Token& Scanner::peek() {
  EnsureTokensInQueue();
  assert(!m_tokens.empty());  // should we be asserting here? I mean, we really
                              // just be checking
                              // if it's empty before peeking.

#if 0
		static Token *pLast = 0;
		if(pLast != &m_tokens.front())
			std::cerr << "peek: " << m_tokens.front() << "\n";
		pLast = &m_tokens.front();
#endif

  return m_tokens.front();
}

// mark
// . Returns the current mark in the stream
Mark Scanner::mark() const { return INPUT.mark(); }

// EnsureTokensInQueue
// . Scan until there's a valid token at the front of the queue,
//   or we're sure the queue is empty.
void Scanner::EnsureTokensInQueue() {
  while (1) {
    if (!m_tokens.empty()) {
      Token& token = m_tokens.front();

      // if this guy's valid, then we're done
      if (token.status == Token::VALID)
        return;

      // here's where we clean up the impossible tokens
      if (token.status == Token::INVALID) {
        m_tokens.pop();
        continue;
      }

      // note: what's left are the unverified tokens
    }

    // no token? maybe we've actually finished
    if (m_endedStream)
      return;

    // no? then scan...
    ScanNextToken();
  }
}

// ScanNextToken
// . The main scanning function; here we branch out and
//   scan whatever the next token should be.
void Scanner::ScanNextToken() {
  if (m_endedStream)
    return;

  if (!m_startedStream)
    return StartStream();

  // get rid of whitespace, etc. (in between tokens it should be irrelevent)
  ScanToNextToken();

  // maybe need to end some blocks
  PopIndentToHere();

  // *****
  // And now branch based on the next few characters!
  // *****

  // end of stream
  if (!INPUT)
    return EndStream();

  if (INPUT.column() == 0 && INPUT.peek() == Keys::Directive)
    return ScanDirective();

  // document token
  if (INPUT.column() == 0 && Exp::DocStart().Matches(INPUT))
    return ScanDocStart();

  if (INPUT.column() == 0 && Exp::DocEnd().Matches(INPUT))
    return ScanDocEnd();

  // flow start/end/entry
  if (INPUT.peek() == Keys::FlowSeqStart || INPUT.peek() == Keys::FlowMapStart)
    return ScanFlowStart();

  if (INPUT.peek() == Keys::FlowSeqEnd || INPUT.peek() == Keys::FlowMapEnd)
    return ScanFlowEnd();

  if (INPUT.peek() == Keys::FlowEntry)
    return ScanFlowEntry();

  // block/map stuff
  if (Exp::BlockEntry().Matches(INPUT))
    return ScanBlockEntry();

  if ((InBlockContext() ? Exp::Key() : Exp::KeyInFlow()).Matches(INPUT))
    return ScanKey();

  if (GetValueRegex().Matches(INPUT))
    return ScanValue();

  // alias/anchor
  if (INPUT.peek() == Keys::Alias || INPUT.peek() == Keys::Anchor)
    return ScanAnchorOrAlias();

  // tag
  if (INPUT.peek() == Keys::Tag)
    return ScanTag();

  // special scalars
  if (InBlockContext() && (INPUT.peek() == Keys::LiteralScalar ||
                           INPUT.peek() == Keys::FoldedScalar))
    return ScanBlockScalar();

  if (INPUT.peek() == '\'' || INPUT.peek() == '\"')
    return ScanQuotedScalar();

  // plain scalars
  if ((InBlockContext() ? Exp::PlainScalar() : Exp::PlainScalarInFlow())
          .Matches(INPUT))
    return ScanPlainScalar();

  // don't know what it is!
  throw ParserException(INPUT.mark(), ErrorMsg::UNKNOWN_TOKEN);
}

// ScanToNextToken
// . Eats input until we reach the next token-like thing.
void Scanner::ScanToNextToken() {
  while (1) {
    // first eat whitespace
    while (INPUT && IsWhitespaceToBeEaten(INPUT.peek())) {
      if (InBlockContext() && Exp::Tab().Matches(INPUT))
        m_simpleKeyAllowed = false;
      INPUT.eat(1);
    }

    // then eat a comment
    if (Exp::Comment().Matches(INPUT)) {
      // eat until line break
      while (INPUT && !Exp::Break().Matches(INPUT))
        INPUT.eat(1);
    }

    // if it's NOT a line break, then we're done!
    if (!Exp::Break().Matches(INPUT))
      break;

    // otherwise, let's eat the line break and keep going
    int n = Exp::Break().Match(INPUT);
    INPUT.eat(n);

    // oh yeah, and let's get rid of that simple key
    InvalidateSimpleKey();

    // new line - we may be able to accept a simple key now
    if (InBlockContext())
      m_simpleKeyAllowed = true;
  }
}

///////////////////////////////////////////////////////////////////////
// Misc. helpers

// IsWhitespaceToBeEaten
// . We can eat whitespace if it's a space or tab
// . Note: originally tabs in block context couldn't be eaten
//         "where a simple key could be allowed
//         (i.e., not at the beginning of a line, or following '-', '?', or
// ':')"
//   I think this is wrong, since tabs can be non-content whitespace; it's just
//   that they can't contribute to indentation, so once you've seen a tab in a
//   line, you can't start a simple key
bool Scanner::IsWhitespaceToBeEaten(char ch) {
  if (ch == ' ')
    return true;

  if (ch == '\t')
    return true;

  return false;
}

// GetValueRegex
// . Get the appropriate regex to check if it's a value token
const RegEx& Scanner::GetValueRegex() const {
  if (InBlockContext())
    return Exp::Value();

  return m_canBeJSONFlow ? Exp::ValueInJSONFlow() : Exp::ValueInFlow();
}

// StartStream
// . Set the initial conditions for starting a stream.
void Scanner::StartStream() {
  m_startedStream = true;
  m_simpleKeyAllowed = true;
  std::auto_ptr<IndentMarker> pIndent(new IndentMarker(-1, IndentMarker::NONE));
  m_indentRefs.push_back(pIndent);
  m_indents.push(&m_indentRefs.back());
}

// EndStream
// . Close out the stream, finish up, etc.
void Scanner::EndStream() {
  // force newline
  if (INPUT.column() > 0)
    INPUT.ResetColumn();

  PopAllIndents();
  PopAllSimpleKeys();

  m_simpleKeyAllowed = false;
  m_endedStream = true;
}

Token* Scanner::PushToken(Token::TYPE type) {
  m_tokens.push(Token(type, INPUT.mark()));
  return &m_tokens.back();
}

Token::TYPE Scanner::GetStartTokenFor(IndentMarker::INDENT_TYPE type) const {
  switch (type) {
    case IndentMarker::SEQ:
      return Token::BLOCK_SEQ_START;
    case IndentMarker::MAP:
      return Token::BLOCK_MAP_START;
    case IndentMarker::NONE:
      assert(false);
      break;
  }
  assert(false);
  throw std::runtime_error("yaml-cpp: internal error, invalid indent type");
}

// PushIndentTo
// . Pushes an indentation onto the stack, and enqueues the
//   proper token (sequence start or mapping start).
// . Returns the indent marker it generates (if any).
Scanner::IndentMarker* Scanner::PushIndentTo(int column,
                                             IndentMarker::INDENT_TYPE type) {
  // are we in flow?
  if (InFlowContext())
    return 0;

  std::auto_ptr<IndentMarker> pIndent(new IndentMarker(column, type));
  IndentMarker& indent = *pIndent;
  const IndentMarker& lastIndent = *m_indents.top();

  // is this actually an indentation?
  if (indent.column < lastIndent.column)
    return 0;
  if (indent.column == lastIndent.column &&
      !(indent.type == IndentMarker::SEQ &&
        lastIndent.type == IndentMarker::MAP))
    return 0;

  // push a start token
  indent.pStartToken = PushToken(GetStartTokenFor(type));

  // and then the indent
  m_indents.push(&indent);
  m_indentRefs.push_back(pIndent);
  return &m_indentRefs.back();
}

// PopIndentToHere
// . Pops indentations off the stack until we reach the current indentation
// level,
//   and enqueues the proper token each time.
// . Then pops all invalid indentations off.
void Scanner::PopIndentToHere() {
  // are we in flow?
  if (InFlowContext())
    return;

  // now pop away
  while (!m_indents.empty()) {
    const IndentMarker& indent = *m_indents.top();
    if (indent.column < INPUT.column())
      break;
    if (indent.column == INPUT.column() && !(indent.type == IndentMarker::SEQ &&
                                             !Exp::BlockEntry().Matches(INPUT)))
      break;

    PopIndent();
  }

  while (!m_indents.empty() && m_indents.top()->status == IndentMarker::INVALID)
    PopIndent();
}

// PopAllIndents
// . Pops all indentations (except for the base empty one) off the stack,
//   and enqueues the proper token each time.
void Scanner::PopAllIndents() {
  // are we in flow?
  if (InFlowContext())
    return;

  // now pop away
  while (!m_indents.empty()) {
    const IndentMarker& indent = *m_indents.top();
    if (indent.type == IndentMarker::NONE)
      break;

    PopIndent();
  }
}

// PopIndent
// . Pops a single indent, pushing the proper token
void Scanner::PopIndent() {
  const IndentMarker& indent = *m_indents.top();
  m_indents.pop();

  if (indent.status != IndentMarker::VALID) {
    InvalidateSimpleKey();
    return;
  }

  if (indent.type == IndentMarker::SEQ)
    m_tokens.push(Token(Token::BLOCK_SEQ_END, INPUT.mark()));
  else if (indent.type == IndentMarker::MAP)
    m_tokens.push(Token(Token::BLOCK_MAP_END, INPUT.mark()));
}

// GetTopIndent
int Scanner::GetTopIndent() const {
  if (m_indents.empty())
    return 0;
  return m_indents.top()->column;
}

// ThrowParserException
// . Throws a ParserException with the current token location
//   (if available).
// . Does not parse any more tokens.
void Scanner::ThrowParserException(const std::string& msg) const {
  Mark mark = Mark::null_mark();
  if (!m_tokens.empty()) {
    const Token& token = m_tokens.front();
    mark = token.mark;
  }
  throw ParserException(mark, msg);
}
}

#include "src/scanscalar.h"

#include <algorithm>

#include "src/exp.h"
#include "src/regeximpl.h"
#include "src/stream.h"
#include "yaml-cpp/exceptions.h"  // IWYU pragma: keep

namespace YAML {
// ScanScalar
// . This is where the scalar magic happens.
//
// . We do the scanning in three phases:
//   1. Scan until newline
//   2. Eat newline
//   3. Scan leading blanks.
//
// . Depending on the parameters given, we store or stop
//   and different places in the above flow.
std::string ScanScalar(Stream& INPUT, ScanScalarParams& params) {
  bool foundNonEmptyLine = false;
  bool pastOpeningBreak = (params.fold == FOLD_FLOW);
  bool emptyLine = false, moreIndented = false;
  int foldedNewlineCount = 0;
  bool foldedNewlineStartedMoreIndented = false;
  std::size_t lastEscapedChar = std::string::npos;
  std::string scalar;
  params.leadingSpaces = false;

  while (INPUT) {
    // ********************************
    // Phase #1: scan until line ending

    std::size_t lastNonWhitespaceChar = scalar.size();
    bool escapedNewline = false;
    while (!params.end.Matches(INPUT) && !Exp::Break().Matches(INPUT)) {
      if (!INPUT)
        break;

      // document indicator?
      if (INPUT.column() == 0 && Exp::DocIndicator().Matches(INPUT)) {
        if (params.onDocIndicator == BREAK)
          break;
        else if (params.onDocIndicator == THROW)
          throw ParserException(INPUT.mark(), ErrorMsg::DOC_IN_SCALAR);
      }

      foundNonEmptyLine = true;
      pastOpeningBreak = true;

      // escaped newline? (only if we're escaping on slash)
      if (params.escape == '\\' && Exp::EscBreak().Matches(INPUT)) {
        // eat escape character and get out (but preserve trailing whitespace!)
        INPUT.get();
        lastNonWhitespaceChar = scalar.size();
        lastEscapedChar = scalar.size();
        escapedNewline = true;
        break;
      }

      // escape this?
      if (INPUT.peek() == params.escape) {
        scalar += Exp::Escape(INPUT);
        lastNonWhitespaceChar = scalar.size();
        lastEscapedChar = scalar.size();
        continue;
      }

      // otherwise, just add the damn character
      char ch = INPUT.get();
      scalar += ch;
      if (ch != ' ' && ch != '\t')
        lastNonWhitespaceChar = scalar.size();
    }

    // eof? if we're looking to eat something, then we throw
    if (!INPUT) {
      if (params.eatEnd)
        throw ParserException(INPUT.mark(), ErrorMsg::EOF_IN_SCALAR);
      break;
    }

    // doc indicator?
    if (params.onDocIndicator == BREAK && INPUT.column() == 0 &&
        Exp::DocIndicator().Matches(INPUT))
      break;

    // are we done via character match?
    int n = params.end.Match(INPUT);
    if (n >= 0) {
      if (params.eatEnd)
        INPUT.eat(n);
      break;
    }

    // do we remove trailing whitespace?
    if (params.fold == FOLD_FLOW)
      scalar.erase(lastNonWhitespaceChar);

    // ********************************
    // Phase #2: eat line ending
    n = Exp::Break().Match(INPUT);
    INPUT.eat(n);

    // ********************************
    // Phase #3: scan initial spaces

    // first the required indentation
    while (INPUT.peek() == ' ' && (INPUT.column() < params.indent ||
                                   (params.detectIndent && !foundNonEmptyLine)))
      INPUT.eat(1);

    // update indent if we're auto-detecting
    if (params.detectIndent && !foundNonEmptyLine)
      params.indent = std::max(params.indent, INPUT.column());

    // and then the rest of the whitespace
    while (Exp::Blank().Matches(INPUT)) {
      // we check for tabs that masquerade as indentation
      if (INPUT.peek() == '\t' && INPUT.column() < params.indent &&
          params.onTabInIndentation == THROW)
        throw ParserException(INPUT.mark(), ErrorMsg::TAB_IN_INDENTATION);

      if (!params.eatLeadingWhitespace)
        break;

      INPUT.eat(1);
    }

    // was this an empty line?
    bool nextEmptyLine = Exp::Break().Matches(INPUT);
    bool nextMoreIndented = Exp::Blank().Matches(INPUT);
    if (params.fold == FOLD_BLOCK && foldedNewlineCount == 0 && nextEmptyLine)
      foldedNewlineStartedMoreIndented = moreIndented;

    // for block scalars, we always start with a newline, so we should ignore it
    // (not fold or keep)
    if (pastOpeningBreak) {
      switch (params.fold) {
        case DONT_FOLD:
          scalar += "\n";
          break;
        case FOLD_BLOCK:
          if (!emptyLine && !nextEmptyLine && !moreIndented &&
              !nextMoreIndented && INPUT.column() >= params.indent)
            scalar += " ";
          else if (nextEmptyLine)
            foldedNewlineCount++;
          else
            scalar += "\n";

          if (!nextEmptyLine && foldedNewlineCount > 0) {
            scalar += std::string(foldedNewlineCount - 1, '\n');
            if (foldedNewlineStartedMoreIndented ||
                nextMoreIndented | !foundNonEmptyLine)
              scalar += "\n";
            foldedNewlineCount = 0;
          }
          break;
        case FOLD_FLOW:
          if (nextEmptyLine)
            scalar += "\n";
          else if (!emptyLine && !nextEmptyLine && !escapedNewline)
            scalar += " ";
          break;
      }
    }

    emptyLine = nextEmptyLine;
    moreIndented = nextMoreIndented;
    pastOpeningBreak = true;

    // are we done via indentation?
    if (!emptyLine && INPUT.column() < params.indent) {
      params.leadingSpaces = true;
      break;
    }
  }

  // post-processing
  if (params.trimTrailingSpaces) {
    std::size_t pos = scalar.find_last_not_of(' ');
    if (lastEscapedChar != std::string::npos) {
      if (pos < lastEscapedChar || pos == std::string::npos)
        pos = lastEscapedChar;
    }
    if (pos < scalar.size())
      scalar.erase(pos + 1);
  }

  switch (params.chomp) {
    case CLIP: {
      std::size_t pos = scalar.find_last_not_of('\n');
      if (lastEscapedChar != std::string::npos) {
        if (pos < lastEscapedChar || pos == std::string::npos)
          pos = lastEscapedChar;
      }
      if (pos == std::string::npos)
        scalar.erase();
      else if (pos + 1 < scalar.size())
        scalar.erase(pos + 2);
    } break;
    case STRIP: {
      std::size_t pos = scalar.find_last_not_of('\n');
      if (lastEscapedChar != std::string::npos) {
        if (pos < lastEscapedChar || pos == std::string::npos)
          pos = lastEscapedChar;
      }
      if (pos == std::string::npos)
        scalar.erase();
      else if (pos < scalar.size())
        scalar.erase(pos + 1);
    } break;
    default:
      break;
  }

  return scalar;
}
}

#include "src/exp.h"
#include "src/regex_yaml.h"
#include "src/regeximpl.h"
#include "src/stream.h"
#include "yaml-cpp/exceptions.h"  // IWYU pragma: keep
#include "yaml-cpp/mark.h"

namespace YAML {
const std::string ScanVerbatimTag(Stream& INPUT) {
  std::string tag;

  // eat the start character
  INPUT.get();

  while (INPUT) {
    if (INPUT.peek() == Keys::VerbatimTagEnd) {
      // eat the end character
      INPUT.get();
      return tag;
    }

    int n = Exp::URI().Match(INPUT);
    if (n <= 0)
      break;

    tag += INPUT.get(n);
  }

  throw ParserException(INPUT.mark(), ErrorMsg::END_OF_VERBATIM_TAG);
}

const std::string ScanTagHandle(Stream& INPUT, bool& canBeHandle) {
  std::string tag;
  canBeHandle = true;
  Mark firstNonWordChar;

  while (INPUT) {
    if (INPUT.peek() == Keys::Tag) {
      if (!canBeHandle)
        throw ParserException(firstNonWordChar, ErrorMsg::CHAR_IN_TAG_HANDLE);
      break;
    }

    int n = 0;
    if (canBeHandle) {
      n = Exp::Word().Match(INPUT);
      if (n <= 0) {
        canBeHandle = false;
        firstNonWordChar = INPUT.mark();
      }
    }

    if (!canBeHandle)
      n = Exp::Tag().Match(INPUT);

    if (n <= 0)
      break;

    tag += INPUT.get(n);
  }

  return tag;
}

const std::string ScanTagSuffix(Stream& INPUT) {
  std::string tag;

  while (INPUT) {
    int n = Exp::Tag().Match(INPUT);
    if (n <= 0)
      break;

    tag += INPUT.get(n);
  }

  if (tag.empty())
    throw ParserException(INPUT.mark(), ErrorMsg::TAG_WITH_NO_SUFFIX);

  return tag;
}
}