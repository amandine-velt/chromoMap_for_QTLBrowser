#!/usr/bin/bash

cd /home/avelt/data2/2024_chromoMap_for_QTLBrowser_Examples
module load python/3.8
python -m http.server 8000

# Go to : 
# http://0.0.0.0:8000/