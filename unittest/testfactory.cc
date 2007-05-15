/* vim: set sw=4 sts=4 et foldmethod=syntax : */
#include "testfactory.hh"


BaseTest * TestFactory::NewTest (const string & description)
{
    if (description == "BlasTest")
        return new BlasTest;
    if (description == "LiblaTest")
        return new LiblaTest;
    return NULL;
}

