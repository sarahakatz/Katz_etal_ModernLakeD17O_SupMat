text
{% include 3dscatter.html %}
htmltools::includeHTML("3dscatter.html")

title: "test"
output:
  html_document:
    includes:
      before_body: 3dscatter.html
