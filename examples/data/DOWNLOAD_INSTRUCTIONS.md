# Download Instructions for Licensed Datasets

These datasets require a free account or license and cannot be
redistributed. Download them manually and place in this directory.

## Amazon India Product Prices
File needed: amazon_prices.csv (columns: name, actual_price_inr)
Source: https://www.kaggle.com/datasets/karkavelrajaj/amazon-sales-dataset
License: CC0 Public Domain
Steps:
  1. Log in to Kaggle
  2. Download amazon.csv
  3. Run: python examples/scripts/prepare_amazon.py amazon.csv

## YouTube Brazil Trending
File needed: youtube_likes.csv (columns: video_id, likes)
Source: https://www.kaggle.com/datasets/rsrishav/youtube-trending-video-dataset
License: CC0 Public Domain
Steps:
  1. Download BR_Trending.csv from Kaggle
  2. Run: python examples/scripts/prepare_youtube.py BR_Trending.csv

## Intel Stock (INTC)
File needed: intel_stock.csv (columns: date, close, volume)
Source: Yahoo Finance (free) or:
        https://www.kaggle.com/datasets (search INTC)
Steps:
  1. Download INTC.csv
  2. Copy to: examples/data/intel_stock.csv

## AnAge Body Masses
File needed: anage_bodymass.csv
Source: https://genomics.senescence.info/species/
License: CC-BY 3.0
Steps:
  1. Download anage_data.zip from the AnAge website
  2. Extract anage_data.txt
  3. Run: python examples/scripts/prepare_anage.py anage_data.txt

## HYG Stellar Catalog
File needed: hyg_luminosities.csv
Source: https://github.com/astronexus/HYG-Database (CC-BY-SA 2.5)
Steps:
  1. Download hygdata_v42.csv.gz from the repository
  2. Run: python examples/scripts/prepare_hyg.py hygdata_v42.csv.gz
