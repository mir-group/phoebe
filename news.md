---
title: News
layout: template
filename: news
---

## News

The code is soon to be released, stay tuned!

<ul>
  {% for post in site.posts %}
    <li>
      <a href="{{ post.url }}">{{ post.title }}</a>
    </li>
  {% endfor %}
</ul>
