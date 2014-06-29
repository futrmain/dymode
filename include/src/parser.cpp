#include <cstdio>
#include <sstream>

#include "directives.h"  // IWYU pragma: keep
#include "scanner.h"     // IWYU pragma: keep
#include "singledocparser.h"
#include "token.h"
#include "yaml-cpp/exceptions.h"  // IWYU pragma: keep
#include "yaml-cpp/parser.h"

namespace YAML {
class EventHandler;

Parser::Parser() {}

Parser::Parser(std::istream& in) { Load(in); }

Parser::~Parser() {}

Parser::operator bool() const {
  return m_pScanner.get() && !m_pScanner->empty();
}

void Parser::Load(std::istream& in) {
  m_pScanner.reset(new Scanner(in));
  m_pDirectives.reset(new Directives);
}

// HandleNextDocument
// . Handles the next document
// . Throws a ParserException on error.
// . Returns false if there are no more documents
bool Parser::HandleNextDocument(EventHandler& eventHandler) {
  if (!m_pScanner.get())
    return false;

  ParseDirectives();
  if (m_pScanner->empty())
    return false;

  SingleDocParser sdp(*m_pScanner, *m_pDirectives);
  sdp.HandleDocument(eventHandler);
  return true;
}

// ParseDirectives
// . Reads any directives that are next in the queue.
void Parser::ParseDirectives() {
  bool readDirective = false;

  while (1) {
    if (m_pScanner->empty())
      break;

    Token& token = m_pScanner->peek();
    if (token.type != Token::DIRECTIVE)
      break;

    // we keep the directives from the last document if none are specified;
    // but if any directives are specific, then we reset them
    if (!readDirective)
      m_pDirectives.reset(new Directives);

    readDirective = true;
    HandleDirective(token);
    m_pScanner->pop();
  }
}

void Parser::HandleDirective(const Token& token) {
  if (token.value == "YAML")
    HandleYamlDirective(token);
  else if (token.value == "TAG")
    HandleTagDirective(token);
}

// HandleYamlDirective
// . Should be of the form 'major.minor' (like a version number)
void Parser::HandleYamlDirective(const Token& token) {
  if (token.params.size() != 1)
    throw ParserException(token.mark, ErrorMsg::YAML_DIRECTIVE_ARGS);

  if (!m_pDirectives->version.isDefault)
    throw ParserException(token.mark, ErrorMsg::REPEATED_YAML_DIRECTIVE);

  std::stringstream str(token.params[0]);
  str >> m_pDirectives->version.major;
  str.get();
  str >> m_pDirectives->version.minor;
  if (!str || str.peek() != EOF)
    throw ParserException(
        token.mark, std::string(ErrorMsg::YAML_VERSION) + token.params[0]);

  if (m_pDirectives->version.major > 1)
    throw ParserException(token.mark, ErrorMsg::YAML_MAJOR_VERSION);

  m_pDirectives->version.isDefault = false;
  // TODO: warning on major == 1, minor > 2?
}

// HandleTagDirective
// . Should be of the form 'handle prefix', where 'handle' is converted to
// 'prefix' in the file.
void Parser::HandleTagDirective(const Token& token) {
  if (token.params.size() != 2)
    throw ParserException(token.mark, ErrorMsg::TAG_DIRECTIVE_ARGS);

  const std::string& handle = token.params[0];
  const std::string& prefix = token.params[1];
  if (m_pDirectives->tags.find(handle) != m_pDirectives->tags.end())
    throw ParserException(token.mark, ErrorMsg::REPEATED_TAG_DIRECTIVE);

  m_pDirectives->tags[handle] = prefix;
}

void Parser::PrintTokens(std::ostream& out) {
  if (!m_pScanner.get())
    return;

  while (1) {
    if (m_pScanner->empty())
      break;

    out << m_pScanner->peek() << "\n";
    m_pScanner->pop();
  }
}
}

#include "src/regex_yaml.h"

namespace YAML {
// constructors
RegEx::RegEx() : m_op(REGEX_EMPTY) {}

RegEx::RegEx(REGEX_OP op) : m_op(op) {}

RegEx::RegEx(char ch) : m_op(REGEX_MATCH), m_a(ch) {}

RegEx::RegEx(char a, char z) : m_op(REGEX_RANGE), m_a(a), m_z(z) {}

RegEx::RegEx(const std::string& str, REGEX_OP op) : m_op(op) {
  for (std::size_t i = 0; i < str.size(); i++)
    m_params.push_back(RegEx(str[i]));
}

// combination constructors
RegEx operator!(const RegEx & ex) {
  RegEx ret(REGEX_NOT);
  ret.m_params.push_back(ex);
  return ret;
}

RegEx operator||(const RegEx& ex1, const RegEx& ex2) {
  RegEx ret(REGEX_OR);
  ret.m_params.push_back(ex1);
  ret.m_params.push_back(ex2);
  return ret;
}

RegEx operator&&(const RegEx& ex1, const RegEx& ex2) {
  RegEx ret(REGEX_AND);
  ret.m_params.push_back(ex1);
  ret.m_params.push_back(ex2);
  return ret;
}

RegEx operator+(const RegEx& ex1, const RegEx& ex2) {
  RegEx ret(REGEX_SEQ);
  ret.m_params.push_back(ex1);
  ret.m_params.push_back(ex2);
  return ret;
}
}

#include <cassert>
#include <memory>

#include "src/exp.h"
#include "src/scanner.h"
#include "src/token.h"
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