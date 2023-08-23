#ifndef EXCEPTION_HPP
#define EXCEPTION_HPP

// Include header files
#include "TypeTraits.hpp"

// Include libraries
#include <iostream>
#include <exception>

/**
 * Class for the construction of a generic exception, to which it is possible to pass a string content.
 * 
 * It is publicly derived by the standard exception.
*/

namespace TVSFCM{
using T = TypeTraits;

class MyException: public std::exception{
public:
    /**
     * Constructor 
     * @param message_ Content of the exception
    */
    explicit MyException(const T::ExceptionType& message_): 
        message(message_){};

    /**
     * Constructor
     * @param message_ Content of the exception
    */
    explicit MyException(const char* message_):
        message(message_){};

    /**
     * what() method for returning the content of the exception
    */
    virtual const char* what() const noexcept {
        return message.c_str();
    };

    /**
     * Destructor
    */
    virtual ~ MyException() = default;

protected:
    T::ExceptionType message;           //! Content of the exception
};

} // end namespace


#endif // EXCEPTION_HPP