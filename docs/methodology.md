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