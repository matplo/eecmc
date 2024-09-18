#!/bin/bash

python -m pip freeze | sed 's/==.*//' > requirements.txt
