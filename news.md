---
title: News
layout: template
filename: news
---

# News

<div class="row">

  {% for post in site.posts %}
  <div class="column">
    <div style="padding:0.25em; border-radius:1rem; margin:0.75em; border:solid 0.25em #d9d9d9;">
      <div style="min-height:1.5em; width:100%; display:block;">

        <div style="text-align:left; width:70%; height:100%; float:left; padding-left:1em;">
        <a href="{{site.baseurl}}/{{ post.url }}"> {{ post.title }}</a>
        </div>
        <div style="text-align:right; width:20%; height:100%; float:right; padding-right:1em;">
        by {{post.author}}
        </div>
      </div>
    </div>
  </div>
  {% endfor %}

</div>

