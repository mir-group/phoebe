@page Website Website

The website is hosted by Github page and is powered by [Jekyll](https://jekyllrb.com/).

In short, with Jekyll we can use Markdown language to write a web page, and we can also add HTML&CSS instructions for better customization of the web page.
Note that it's a static website, i.e. it doesn't have any database attached to it; however, that shouldn't be a problem for our purposes.
Instructions for Jekyll can be found on their [website](https://jekyllrb.com/).


Note for Jekyll posts: Jekyll 3.*, currently used by github at the time of writing, requires a baseurl to link posts correctly (i.e. the string "/phoebe/" after github.io).
For this reason, we added the baseurl in the links to posts in file `news.md`.
After github pages will update Jekyll, these links will likely break, as Jekyll 4.* doesn't require this baseurl.
To fix it, simply remove the baseurl from `news.md`.



@section posts How to add a post

First, we switch to the github page branch of the phoebe repository:
~~~~~~~~~~~~~~~~~~~{.c}
git checkout gh-pages
~~~~~~~~~~~~~~~~~~~
You should see that most files disappeared, as this section contains only the website.
Next go to this folder:
~~~~~~~~~~~~~~~~~~~{.c}
cd _posts
~~~~~~~~~~~~~~~~~~~
Let's create a new file, in the format `"YYY-MM-DD-NAME.md"`, named using the current date and a name for the post.
Next, we edit this file's content, which may look like this:
~~~~~~~~~~~~~~~~~~~{.c}
---
layout: post
title:  "Welcome to Phoebe News!"
author: Phoebe Team
permalink: post_1
---

This is the post content
~~~~~~~~~~~~~~~~~~~

The first part is the so-called Jekyll front-matter.
The variable `layout` should not be modified, as it describes the HTML style to display the post;
`title` is the title of the Post as it will be displayed on the web;
`author` to identify who's writing this post.
<blockquote>
The variable `permalink` is a name for the url and should be unique, make sure you increment this number compared to the previous post!
</blockquote>

Finally, you can write the post content using standard markup language.




@section know Things to know
The files `*.md` are the pages that you can find in the website, with `index.md` being the home page.
Posts must go in the folder `_posts` as required by Jekyll.
The folder `stylesheets` contain the HTML styles used by the Jekyll base template (things like the font, the header section, etc).

The folder `_layouts` contains a description of the common features of website pages.
In particular, `template.html` contains the overall layout of the HTML page, and has been customized adding a top-bar navigation menu.
The file `post.html` inherits from `template.html`, but adds a few characteristics typical of a Post (title, date, author, etc.).





@section newpage How to add a new page

If you want to add a new page, create a new `NEWNAME.md` file.
The file should look like this:
~~~~~~~~~~~~~~~{.c}
---
title: TITLE
layout: template
filename: NEWNAME
---

Page content
~~~~~~~~~~~~~~~
You should replace `TITLE` and `NEWNAME` with suitable new names.
The page content can be written in Markdown language, or you can also write directly in HTML (as you can see e.g. in `team.md`.
Note that you can also add CSS specific to just this page.

You still need to provide links to this page!
If you want to add it to the navigation bar, open the file `_layouts/template.html` and look for the tag:
~~~~~~~~~~~{.c}
<ul class="menu">
~~~~~~~~~~~
This tag contains a list of links to the various pages in the menu.
For example:
~~~~~~~~~~~~~~{.c}
<li class="item"><a style="text-decoration: none" href="team">Team</a></li>
~~~~~~~~~~~~~~
Here, you may copy this statement and modify it with the proper `href` (the path to the file `team.md` relative to the base repository directory) and the page name (here `Team`).






@section website How to run the website locally

If you plan on doing more extensive changes, you'll probably want to run the website locally, rather than pushing to github and wait for the website to be updated.
To this aim, follow these instructions that have been tested on Ubuntu 20.04.

We need to have `ruby` on our computer:
~~~~~~~~~~~~~~~{.c}
sudo apt install ruby-full
~~~~~~~~~~~~~~~

Next, we install two applications: `Jekyll` (the site generator) and `bundler` (which allows us to launch the website / ruby application):
~~~~~~~~~~~~~~~{.c}
gem install bundler jekyll
~~~~~~~~~~~~~~~

Now, let's cd into the /phoebe repository base folder, where you should see a file called `Gemfile`.
This file contains a list of Ruby applications that we need to run the website.
We install these dependencies by typing:
~~~~~~~~~~~~~~~{.c}
bundle install
~~~~~~~~~~~~~~~

Finally, we launch the website on localhost by typing 
~~~~~~~~~~~~~~~{.c}
bundle exec jekyll serve
~~~~~~~~~~~~~~~
Leave this terminal running for as long as you want the website to be available.
The output to terminal should look something like this:
~~~~~~~~~~~~~~~{.c}
Configuration file: /home/cepe/git/phoebe-pages/_config.yml
            Source: /home/cepe/git/phoebe-pages
       Destination: /home/cepe/git/phoebe-pages/_site
 Incremental build: disabled. Enable with --incremental
      Generating... 
                    done in 0.063 seconds.
 Auto-regeneration: enabled for '/home/cepe/git/phoebe-pages'
    Server address: http://127.0.0.1:4000
  Server running... press ctrl-c to stop.
~~~~~~~~~~~~~~~

Now, leaving this command running in terminal, copy the link at the `Server address` that you can find in the command terminal output, paste it in a browser and you can see the website.


