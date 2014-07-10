#pragma once
#ifndef PEIGEN_MPI_TYPE_H
#define PEIGEN_MPI_TYPE_H


namespace peigen
{
	
	template <typename c_type>
	MPI::Datatype MPI_type()
	{	// No MPI_type, default
		return 0;
	}

	template <>
	MPI::Datatype MPI_type<char>()
	{
		return MPI::CHAR;
	}

	template <>
	MPI::Datatype MPI_type<int>()
	{
		return MPI::INT;
	}

	template <>
	MPI::Datatype MPI_type<bool>()
	{
		return MPI::BOOL;
	}

	template <>
	MPI::Datatype MPI_type<float>()
	{
		return MPI::FLOAT;
	}

	template <>
	MPI::Datatype MPI_type<double>()
	{	// No MPI_type, default
		return MPI::DOUBLE;
	}

	template <>
	MPI::Datatype MPI_type<complex<float>>()
	{
		return MPI::COMPLEX;
	}

	template <>
	MPI::Datatype MPI_type<complex<double>>()
	{
		return MPI::DOUBLE_COMPLEX;
	}

	//MPI::CHAR	 char	 char
	//	MPI::WCHAR	 wchar_t	 wchar_t
	//	MPI::SHORT	 signed short	 signed short
	//	MPI::INT	 signed int	 signed int
	//	MPI::LONG	 signed long	 signed long
	//	MPI::SIGNED_CHAR	 signed char	 signed char
	//	MPI::UNSIGNED_CHAR	 unsigned char	 unsigned char
	//	MPI::UNSIGNED_SHORT	 unsigned short	 unsigned short
	//	MPI::UNSIGNED	 unsigned int	 unsigned int
	//	MPI::UNSIGNED_LONG	 unsigned long	 unsigned long int
	//	MPI::FLOAT	 float	 float
	//	MPI::DOUBLE	 double	 double
	//	MPI::LONG_DOUBLE	 long double	 long double
	//	MPI::BOOL		 bool
	//	MPI::COMPLEX		 Complex<float>
	//	MPI::DOUBLE_COMPLEX		 Complex<double>
	//	MPI::LONG_DOUBLE_COMPLEX		 Complex<long double>
	//	MPI::BYTE
	//	MPI::PACKED

}	// end namespace peigen

#endif // PEIGEN_MPI_TYPE_H