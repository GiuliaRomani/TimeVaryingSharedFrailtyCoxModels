#ifndef EXCEPTION_HPP
#define EXCEPTION_HPP

// Include libraries
#include <iostream>
#include <exception>

// Include header files
#include "TypeTraits.hpp"

using T = TypeTraits;

/**
 * Class for the construction of a generic exception, for which it is possible to pass a string content.
 * 
 * It is publicly derived by the standard exception
*/

class MyException: public std::exception{
public:
    /**
     * Constructor 
     * @param message_ Content of the exception
    */
    MyException(const T::ExceptionType& message_): 
        message(message_){};

    /**
     * Constructor
     * @param message_ Content of the exception
    */
    MyException(const char* message_):
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
    virtual ~ MyException() noexcept = default;

protected:
    T::ExceptionType message;           //! Content of the exception
};


#endif // EXCEPTION_HPP