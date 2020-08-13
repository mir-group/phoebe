@page Tests Unit testing
Unit testing helps us  making sure that code functionalities work correctly, giving us more confidence that the code will work even after code modifications, refactoring, additions, etc.
When you see a part of the code that lacks unit-testing, you are encouraged to add tests. When you fix a bug, you are encouraged to add a test for making sure that the bug will be checked for in the future.

We are using the GoogleTest library for unit testing.


@section run How to run tests
In the `build` directory, simply type `make runTests` to generate an executable `runTests`.
Execute `runTests` to launch the available unit tests.


@section add How to add a new test
To add a new test, simply:
* create a `*.cpp` file in the `test` directory (can create subfolders too)
* in this new file, write a TEST object following the Googletest syntax.
As an example, consider this test to verify that the square root is correctly computed:
~~~~~~~~~~~~~~~~{.c}
TEST (SquareRootTest, ZeroAndNegativeNos) { 
  ASSERT_EQ (0.0, squareRoot (0.0));
  ASSERT_NEAR(2.,squareRoot(4.),1.0e-7);
}
~~~~~~~~~~~~~~~~
Here, `TEST` is the class that manages the unit test.
The first argument `SquareRootTest` is the name of the test suite, and the second argument `ZeroAndNegativeNos` is the test name. The function `squareRoot` is (supposed to be) defined in the phoebe source code. We are using the Googletest macro `ASSERT_EQ` and `ASSERT_NEAR` to verify the expected outcome of the function.

Note that Googletest will automatically look for any `TEST` object added to the `test` folder and add it to the `runTests` executable.

A coincise guide to Googletest, describing its primary functionalities, is available at its [GitHub page](https://github.com/google/googletest/blob/master/googletest/docs/primer.md).