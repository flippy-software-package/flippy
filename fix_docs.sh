#!/bin/bash

npx embedme ./README.md || exit 1
bash make_doxygen_html.sh


