@page Start Getting started

@tableofcontents

@section Download Download
This is the first paragraph. Separate paragraphs with empty lines.

This is the second paragraph.

@section Installation Installation

@subsection Linux Linux

@subsection MacOS MacOS

This is List:
- First item.
Bla bla bla bla bla bla bla bla bla bla bla

- Second item (use tab for nested items)
  + nested item 1
  + nested item 2
  + nested item 3

Numbered list:

1. First item
2. Second item

This is normal paragraph.

  This is code block.

Now we will use horizontal ruler:
- - - - - - - - - - -

This will give you italic text: *it is important*.

This will give you bold text: **it is very important**

To introduce a span of code, use backticks: for example, `printf()` function.

Now we can add table with certain alignment:

|Right|Center|Left|
|----:|:----:|:---|
|1    |2     |3   |

We can also add code block in C style:
~~~~~~~~~~~~~~~~~~~~~~~~~~~{.c}
void Animal::set_maturity_age( age_type r )
{
    maturity_age_ = r;
}
~~~~~~~~~~~~~~~~~~~~~~~~~~~

In case of problems, [Google] (http://google.com) will help. More details see in section [Download] (@ref Download).

We can also add a figure:
![Fig.1: Logo] (../Logo.png "Logo")

In case of questions, please contact: <mailto:mail@mail.com> or check <http://www.google.com>

Now equation in text: \f$E=mc^2\f$;

On separate line:
\f[
a^2+b^2=c^2
\f]
