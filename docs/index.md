
👉 with **MyST’s `list` or `card`-style navigation**.

---

# 🏠 Correct `index.md` (working version)

```markdown
---
title: "QIIME2 Under the Hood"
description: "Deep dives into QIIME2 internals: demux, DADA2, and debugging."
---

# 🔬 QIIME2 Under the Hood

Most tutorials show you *what commands to run*.  
This blog focuses on:

👉 **What actually happens underneath.**

---

## 📚 Articles

- [Demultiplexing deep dive](demux.md)
- [DADA2 under the hood](dada2.md)
- [Debugging QIIME2](debugging.md)

---

## 🧪 Example: Debugging QIIME2

```bash
python -m pdb -m q2cli dada2 denoise-single ...