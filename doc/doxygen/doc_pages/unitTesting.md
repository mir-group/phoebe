@page Tests Testing

@section unittest Unit testing

Unit testing helps us  making sure that code functionalities work correctly, giving us more confidence that the code will work even after code modifications, refactoring, additions, etc.
When you see a part of the code that lacks unit-testing, you are encouraged to add tests. When you fix a bug, you are encouraged to add a test for making sure that the bug will be checked for in the future.

We are using the GoogleTest library for unit testing.


@subsection utrun How to run tests
In the `build` directory, simply type
~~~~~~~~~~~~~~~{.c}
make -j4 runTests
~~~~~~~~~~~~~~~
to generate an executable `runTests`.
Execute
~~~~~~~~~~~~~~{.c}
./runTests
~~~~~~~~~~~~~~
to launch the available unit tests.

When creating a new test, you may want to only run the test that matters to you.
For example, type
~~~~~~~~~~~~{.c}
./runTests --gtest_filter=Points*
~~~~~~~~~~~~
to execute all tests in the group Points* (the first argument of the TEST() object initialization, see below).

@subsection utadd How to add a new test
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





@section citest Continuous integration
We are using Github Actions for the continuous integration of the repository.
Everytime a commit is created, the code will do:

<ol>
<li> Download the phoebe repository
<li> Compile the code (and also runTests)
<li> Download data files needed to run the examples.
     Currently, these data files are stored in the phoebe-data github repository
<li> Run the unit tests
<li> Run the Silicon phonon examples and check the output json against some reference json files. The reference json files have been created by phoebe as well (are not experimental data).
<li> Same also for the examples of Silicon electronic transport with Wannier interpolation and EPA.
</ol>

Note also that the CI is tested on a Ubuntu virtual machine, with and without MPI, with and without OpenMP parallelization levels.
We cannot test phoebe on a GPU with the CI (without paying).

The configuration file of Github Actions is stored in
~~~~~~~~{.c}
./.github/workflows/buildandtest.yml
~~~~~~~~
Modify it with caution, if necessary!
This may be needed when adding new features and new CI example-tests.
