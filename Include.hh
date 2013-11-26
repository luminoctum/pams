#ifndef INCLUDE 
#define INCLUDE

/* define error messages */
#define ALARM(message)  std::cout << message << std::endl 
#define ASSERT_FILE_NOT_FOUND(AFILE) \
        std::cerr << "***** ERROR *****" << std::endl; \
        std::cerr << "can not open file \"" << AFILE << "\" ..."<< std::endl; \
        std::cerr << "*****************" << std::endl; \
        assert(0) 
#define ASSERT_VARIABLE_OUT_OF_RANGE(AVAR) \
        std::cerr << "***** ERROR *****" << std::endl; \
        std::cerr << "variable \"" << AVAR << "\" out of range ..."<< std::endl; \
        std::cerr << "*****************" << std::endl; \
        assert(0) 
#define ASSERT_DIMENSION_MISMATCH(AVAR1, AVAR2) \
        std::cerr << "***** ERROR *****" << std::endl; \
        std::cerr << "dimension mismatch for \"" << AVAR1 << "\", \"" << AVAR2 << std::endl; \
        std::cerr << "*****************" << std::endl; \
        assert(0) 
#define ASSERT_NOT_SUPPORTED \
        std::cerr << "***** ERROR *****" << std::endl; \
        std::cerr << "method not supported!" << std::endl; \
        std::cerr << "*****************" << std::endl; \
        assert(0) 
#define BREAKPOINT assert(0)

/* define constants */
#define ZERO1(n) Eigen::ArrayXf::Zero(n)
#define ZERO2(n,m) Eigen::ArrayXXf::Zero(n,m)

/* define small functions */
#define MIN3(a, b, c) ( a < b ? (a < c ? a : c) : (b < c ? b : c) )
#define MIN2(a, b) ( a < b ? a : b )
#define MAX2(a, b) ( a > b ? a : b )

/* define domain decomposition */
#define NTILES_IN_X 1
#define NTILES_IN_Y 1

/* define interpolation order */
#define O_INTERP    2

#endif
