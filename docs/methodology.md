# Methodology

The content here is compiled from my own deep-dives into practical workflows behind the source code, as well as my interpretations of the algorithmic and biological concepts found in official documentation and published research.

Entry interface of QIIME2 codebase is `q2cli`, it employed package [click](https://click.palletsprojects.com/en/stable/) to parse arguments on command line, then decide which required modules for the task. Each module as a wraper, it would decide processing arguments and call needed APIs through system calls (such as `my_api.R arg1 arg2`). For example, when running denoise with DADA2, `q2_dada2` interface of QIIME2 would run [DADA2 R code](https://github.com/benjjneb/dada2) with accompany arguments.

In order to determine which is main function which actually handles the task, debugger tools is employed:

| Debugger tools |   Programing language |
| :------------- |   :------------------ |
| built-in pdb            |   Python              |
| gdb            |   C/C++               |
| built-in R debugger | R                | 
| bashdb              | Bash      |

Debuger helps to set breakpoint and stop a stack (function), line of code and view local variables in runtime. By tracking at the moment, QIIME2 assign task, the the main processing functions would be revealed. Important breakpoints of processes can be found in [breakpoints](https://github.com/Truongphi20/pbwt/blob/master/docs/static/breakpoints).

When main function is dentified, combing with reading docummentation and article to catch the general concepts, and purpose. Base on the references, dissecting it into clear steps of the process and understanding meaning of these. Recurring the process until fully capturing the whole picture.  

It is frequently QIIME2 call a system-call API to handle the task. In these cases, I copy full command tranfered to API to another script and continue to debug as a different sesson. See bash script in [debugging_commands](https://github.com/Truongphi20/qiime2_blog/tree/main/commands) for example. 

The devcontainer with the base image `quay.io/qiime2/amplicon:2026.1` is utilized as a environment to walkthrough QIIME2 codebase, and deploy this blog. See [Dockerfile](https://github.com/Truongphi20/qiime2_blog/blob/main/.devcontainer/Dockerfile) for more detail. 