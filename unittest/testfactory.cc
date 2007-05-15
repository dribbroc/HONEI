/* vim: set sw=4 sts=4 et foldmethod=syntax : */
#include "testfactory.hh"


BaseTest * TestFactory::NewTest (const string & description)
{
    if (description == "BlasTest")
        return new BlasTest("foo");
    if (description == "LiblaTest")
        return new LiblaTest("foo");
    return NULL;
}

