<p align="center"><img src="http://svgur.com/i/JSJ.svg" width=200 ></p>

STK is the **S**PECC **T**ool**k**it, a collection of scripts and tools that make working with the data produced by [SPECC](https://xgitlab.cels.anl.gov/gplynch/specc) easier.

Usage
------

In order to use these functions, first clone this repository using 
```
git clone https://github.com/gplynch619/stk.git
```

If you only plan on using these functions occasionally, add the directory *containing* `stk/` to your path in your import statements:

```
import sys
sys.path.insert(0, "/path/containing/stk")
``` 

If you want to permanently add them to your path, add the containing directory to your `$PYTHONPATH$` environment variable. In `~/.profile` add:
```
export PYTHONPATH=/path/containing/stk:$PYTHONPATH
```

Then you should be able to simply `import stk` as needed.

Prerequisites
------

These functions currently only require an environment with `numpy` and `scipy`. 
