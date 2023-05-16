**This is a development doc that explains code nuances**

Overview
---
gtdb-backend provides REST API to query databases. It supports tree export, search, autocomplete features required by frontend.

It is a python-flask application, with all end points marked by `@app.route` decorator.

**Why do we need @crossdomain decorator?**

The reason is that frontend might be served in a different domain than backend. Although generally, allowing cross origin request for backend API is a bad practice, but in our case, it is not security sensitive, and allows third-party to develop a webpage that uses our data.
