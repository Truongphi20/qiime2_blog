# Methodology

The content presented in this blog is compiled from deep-dives into the practical workflows underlying, alongside my interpretations of the algorithmic and biological concepts found in official documentation and published research.

## QIIME 2 Architecture and Execution
The entry point for the QIIME 2 codebase is `q2cli`, which utilizes the [Click](https://click.palletsprojects.com/en/stable/) package to parse command-line arguments and identify the required modules for a task. Each module acts as a wrapper that processes these arguments and invokes the necessary APIs through system calls (e.g., `my_api.R arg1 arg2`). For example, when performing denoising with DADA2, the `q2-dada2` interface executes the [DADA2 R code](https://github.com/benjjneb/dada2) with the corresponding arguments.

## Identifying Main Functions via Debugging
To determine the primary functions responsible for specific tasks, I employ various debugging tools based on the programming language:

| Debugger Tool | Programming Language |
| :--- | :--- |
| **pdb** (built-in) | Python |
| **gdb** | C/C++ |
| **R debugger** (built-in) | R |
| **bashdb** | Bash |

Debuggers allow for setting breakpoints to pause the call stack or a specific line of code, enabling the inspection of local variables at runtime. By tracing the moment QIIME 2 assigns a task, the core processing functions are revealed. Important breakpoints for these processes are documented in my [breakpoint graphs](https://github.com/Truongphi20/qiime2_blog/blob/main/docs/static/breakpoints).

Sometimes, QIIME 2 utilizes [rachis - former QIIME2 framework](https://github.com/rachis-org/rachis) to embed helper functions from various ecosystem extensions. The function called can be traced by inspecting the `plugin_setup.py` file within each plugin's source code.

## Conceptual Dissection and API Bridges
Once a main function is identified, I combine code analysis with official documentation and research articles to grasp the general purpose and underlying concepts. Based on these references, I dissect the function into clear procedural steps to capture the full picture.

Because QIIME 2 frequently uses system calls to handle tasks, I often encounter "API bridges". In these instances, I copy the full command transferred to the API into a separate script and initiate a new debugging session. Examples of these extracted commands can be found in the [debugging_commands](https://github.com/Truongphi20/qiime2_blog/tree/main/commands) directory.

## Environment and Deployment
A **DevContainer** based on the `quay.io/qiime2/amplicon:2026.1` image is used as the environment for walking through the codebase and deploying this blog. Detailed configurations are available in the [Dockerfile](https://github.com/Truongphi20/qiime2_blog/blob/main/.devcontainer/Dockerfile).