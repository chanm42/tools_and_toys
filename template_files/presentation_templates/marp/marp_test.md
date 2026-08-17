---
marp: true
theme: fzj
math: mathjax
paginate: false
footer: " "
---

<!-- _class: title -->
<!-- _paginate: false -->
# Title
## Sub-title with longer explanation
Author | Date

---

# Body slide - primary part of the presentation
## Body slide subtitle - here we use the body-2col class for automatic 2 column text

<div class="body-2col">

Body slide text block without any bullets
- You can add bullets using - with a space
    - This is a sub-bullet
- This is another main bullet
1. You can make numbered bullet lists using 1.
1. It will automatically number everything using 1.
You can add code stylings as well
```python
def my_func(params):
    x = 1
    y = 2
    return x*y
```

</div>

---

# Example of a body slide with columns using \<div> classes
<div class="columns">

<div class="text">

This is a left side column with some bullet points and equations
- and bullet points
- and a bullet point with equations $H\psi = E\psi$

You can even include much longer sentences and equations $E^2 = P^2c^2 + (mc^2)^2$

</div>

<div class="image">
<figure>

![width:500px height:400px](./img/solar.jpg)
<figcaption>Figure 1: This is a caption for an image using the figure and figcaption containers</figcaption>

</figure>
</div>
</div>

---

# Example of a slide with two images side by side using \<div> elements

<div class="columns">
<div class="image">
<figure>

![width:500px height:400px](./img/solar.jpg)
<figcaption>Figure 2: This image is on the left because it comes first in the div list</figcaption>

</figure>
</div>
<div class="image">
<figure>

![width:500px height:400px](./img/solar.jpg)
<figcaption>Figure 3: This image is on the right because it comes second in the div list</figcaption>

</figure>
</div>
</div>

---

# Example of a slide with two text columns side by side using \<div> elements

<div class="columns">
<div class="text">

This text is in the left side column
- Bullet one
    - Sub bullet one
- Bullet two
    1. Sub bullet list one
    1. Sub bullet list two

</div>
<div class="text red-rectangle">

This text is on the right side of the slide because it is in the second div class

$$\alpha $$

</div>
</div>

---

<!-- _class: transition -->
<!-- _paginate: true -->
<div class="transition-text">

# Here is a transition slide
## You can add sub lines to the transition slide

</div>

<div class="transition-image">

![width=100% height=500px](./img/solar.jpg)

</div>

---

<!-- _class: image -->
# This slide is a large single image with a caption

<div class="image-content">

![./img/solar.jpg](./img/solar.jpg)

</div>
Figure 4: Example caption for figure 4

<script src="./helper.js"></script>
