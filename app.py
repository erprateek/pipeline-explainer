# pipeline_explainer_app/app.py

import streamlit as st
import openai
from pathlib import Path

# === CONFIG ===
openai.api_key = st.secrets["OPENAI_API_KEY"]

# === HELPER FUNCTIONS ===
def split_script_into_blocks(script_text):
    """
    Split script into logical blocks.
    For Snakemake: split by 'rule'
    For bash: split by function or comment blocks
    """
    blocks = []
    lines = script_text.splitlines()
    current_block = []
    for line in lines:
        if line.strip().startswith("rule ") or line.strip().startswith("function ") or line.strip().startswith("#"):
            if current_block:
                blocks.append("\n".join(current_block))
                current_block = []
        current_block.append(line)
    if current_block:
        blocks.append("\n".join(current_block))
    return blocks


def explain_code_block(code_block):
    """Call OpenAI to explain the code block."""
    prompt = f"""
You are a senior bioinformatics scientist.
Explain the following code block to a junior scientist in plain English.
Include relevant tools, inputs/outputs, and what it accomplishes.

Code block:
"""
    prompt += code_block
    response = openai.ChatCompletion.create(
        model="gpt-4-turbo",
        messages=[{"role": "user", "content": prompt}],
        temperature=0.3
    )
    return response.choices[0].message.content.strip()


# === STREAMLIT UI ===
st.set_page_config(page_title="Bioinformatics Pipeline Explainer")
st.title("🧬 Bioinformatics Pipeline Explainer")

uploaded_file = st.file_uploader("Upload a pipeline script (.sh, .smk, .nf)", type=["sh", "smk", "nf", "txt"])

if uploaded_file is not None:
    script_text = uploaded_file.read().decode("utf-8")
    st.subheader("📄 Uploaded Script")
    st.code(script_text, language="bash")

    blocks = split_script_into_blocks(script_text)
    st.subheader("🧠 Explanation")

    for i, block in enumerate(blocks):
        with st.expander(f"Step {i+1}"):
            st.code(block, language="bash")
            explanation = explain_code_block(block)
            st.markdown(explanation)

