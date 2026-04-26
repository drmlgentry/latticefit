import os
os.makedirs("C:/dev/latticefit/.streamlit", exist_ok=True)
cfg = "[theme]\nprimaryColor = \"#4a90d9\"\nbackgroundColor = \"#ffffff\"\nsecondaryBackgroundColor = \"#f4f8ff\"\ntextColor = \"#222222\"\nfont = \"sans serif\"\n"
open("C:/dev/latticefit/.streamlit/config.toml", "w").write(cfg)
print("Done.")
