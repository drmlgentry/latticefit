import re

with open("latticefit_app.py", encoding="utf-8") as f:
    src = f.read()

# Fix 1: Update model name
src = src.replace(
    'model="claude-sonnet-4-20250514"',
    'model="claude-sonnet-4-5"'
)

# Fix 2: Graceful API key handling
old = '''def call_claude(system_prompt, user_message, max_tokens=800):
    if not HAS_ANTHROPIC:
        return "Install anthropic package: pip install anthropic"
    client = anthropic.Anthropic()'''

new = '''def call_claude(system_prompt, user_message, max_tokens=800):
    if not HAS_ANTHROPIC:
        return "Install anthropic package: pip install anthropic"
    import os
    api_key = os.environ.get("ANTHROPIC_API_KEY", "")
    if not api_key:
        return ("AI assistant requires ANTHROPIC_API_KEY environment variable. "
                "Set it with: $env:ANTHROPIC_API_KEY = 'YOUR_API_KEY_HERE'")
    client = anthropic.Anthropic(api_key=api_key)'''

src = src.replace(old, new)

with open("latticefit_app.py", "w", encoding="utf-8") as f:
    f.write(src)
print("Fixed.")
