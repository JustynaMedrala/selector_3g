## Running the code

To compile and run the program, use the following commands:

```bash
root -l -b -q compile.C
root -l -b -q run.C
```

## Changing the input file path

To use a different input file, edit the path in `run.C`:

```cpp
Selector3g s(
    "/data/4/users/jsowa/data/Run10/MC/2025_10_16-23_08_22.ntu.root", // input file
    "output.root", // output file
);
```
## Output

- The main output file is `output.root`, which contains the selected 3 + 1 gamma events from the input file.
- A log file `cut_analysis.log` is generated, containing a summary of all cuts applied.
