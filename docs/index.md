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
```

:::{attention} Disclaim

I am not a member of the QIIME 2 core development team, just a curious bioinformatician wanting to understand what happens under the hood. I’ve written this blog as a reference for my future self and for anyone else sharing this curiosity. 

The content here is compiled from my own deep-dives into workflows behinds the source code, as well as my interpretations of the algorithmic and biological concepts found in official documentation and published research.

:::