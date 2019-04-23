# dmf-chip #

Digital microfluidic chip tools

<!-- vim-markdown-toc GFM -->

* [Command-line API](#command-line-api)
    * [Command `info`](#command-info)
        * [Example](#example)
    * [Command `neighbour detect`](#command-neighbour-detect)
    * [Command `route compute`](#command-route-compute)
        * [Example](#example-1)
* [Command `route list`](#command-route-list)
        * [Example](#example-2)
* [Command `route add`](#command-route-add)
        * [Example](#example-3)
* [Command `route remove`](#command-route-remove)
        * [Example](#example-4)
    * [Install](#install)
    * [License](#license)
    * [Contributors](#contributors)

<!-- vim-markdown-toc -->

-------------------------------------------------------------------------------

# Command-line API

```sh
Usage: dmf-chip [OPTIONS] COMMAND [ARGS]...

Options:
  -v, --verbosity LVL  Either CRITICAL, ERROR, WARNING, INFO or DEBUG
  --help               Show this message and exit.

Commands:
  info       Display information about digital...
  neighbour  Neighbour actions
  route      Route actions
```

## Command `info`

```
Usage: dmf-chip info [OPTIONS] CHIP_FILE

  Display information about digital microfluidics chip design.

Options:
  --output FILENAME
  --help             Show this message and exit.
```

### Example

```sh
$ dmf-chip info chip.svg
Detected Inkscape 0.92+; using 96 pixels per inch
File:                    chip.svg
SHA256:                  c0928207a3cb19841342e4721b8bd6fd14e9f3f4ca5d2df4859caa62baada19e
Pixels per inch:         96
Number of electrodes:    84
Number of channels used: 84
At least one path exists between all electrodes.
```

## Command `neighbour detect`

```sh
Usage: dmf-chip neighbour detect [OPTIONS] CHIP_FILE [OUTPUT]

  Detect adjacent electrodes and add "Connections" layer  to output SVG file
  containing lines; each corresponding to a cached connection between two
  electrodes.

Options:
  --distance-threshold TEXT  Default: 0.175 millimeter
  --help                     Show this message and exit.
```

## Command `route compute`

```sh
Usage: dmf-chip route compute [OPTIONS] CHIP_FILE [OUTPUT]

  Find low-cost test tour that passes through each electrode **at least**
  once and add test route to `<dmf:`ChipDesign><dmf:TestRoutes>`.

Options:
  --start-id TEXT          Start electrode id
  --start-channel INTEGER  Start channel number
  --seed INTEGER           Random seed
  --route-id TEXT          Test route id
  --help                   Show this message and exit.
```

### Example

Compute low-cost test tour for an input file named `input-chip.svg`:

```sh
$ dmf_chip route compute input-chip.svg --route-id default --start-id reservoir-A0 output.svg
Detected Inkscape 0.92+; using 96 pixels per inch
Found 3-opt with gain 8.0
Found 3-opt with gain 2.0
Found 3-opt with gain 2.0
Found 3-opt with gain 2.0
Found 2-opt with gain 2.0
Found 4-opt with gain 4.0
Found 3-opt with gain 2.0
Found 4-opt with gain 2.0
Found 2-opt with gain 2.0
Found 5-opt with gain 2.0
start: 8, cost: 108.0
Found 3-opt with gain 8.0
Found 4-opt with gain 8.0
Found 3-opt with gain 2.0
Found 12-opt with gain 2.0
Found 2-opt with gain 2.0
Found 2-opt with gain 2.0
Found 4-opt with gain 2.0
Found 3-opt with gain 2.0
start: 18, cost: 110.0
Found 5-opt with gain 2.0
Found 3-opt with gain 2.0
Found 3-opt with gain 2.0
Found 4-opt with gain 4.0
Found 3-opt with gain 2.0
Found 3-opt with gain 2.0
start: 57, cost: 112.0
Found 9-opt with gain 2.0
Found 5-opt with gain 4.0
Found 3-opt with gain 2.0
Found 3-opt with gain 2.0
Found 4-opt with gain 2.0
Found 3-opt with gain 2.0
Found 2-opt with gain 2.0
start: 66, cost: 108.0
Found 3-opt with gain 8.0
Found 3-opt with gain 2.0
Found 4-opt with gain 12.0
Found 3-opt with gain 2.0
start: 111, cost: 110.0
found element: `dmf:ChipDesign`
found element: `dmf:TestRoutes`
Add new element: `dmf:TestRoute[@version="0.1.0"][@id="default"]`
Added 84 waypoints.
Wrote test route to: `output.svg`
```

# Command `route list`

```sh
Usage: dmf-chip route list [OPTIONS] CHIP_FILE

  Display routes

Options:
  --help  Show this message and exit.
```

### Example

List test routes in `chip.svg`:

```sh
$ dmf-chip route list chip.svg
Detected Inkscape 0.92+; using 96 pixels per inch
{'version': '0.1.0', 'id': 'default', 'waypoints': ['110', '93', '85', '70', '63', '62', '118', '1', '57', '56', '49', '34', '26', '9', '0', '119']}
{'version': '0.1.0', 'id': 'default2', 'waypoints': ['0', '10']}
```

# Command `route add`

```sh
Usage: dmf-chip route add [OPTIONS] CHIP_FILE [OUTPUT]

  Add a test route

Options:
  --id TEXT               Route id  [required]
  -w, --waypoint INTEGER  Channel number. Multiple channel numbers may be
                          specified, e.g., `-w 10 -w 23 ...`
  --help                  Show this message and exit.
```

### Example

Add test route with `id="default"` to `chip.svg`:

```sh
$ dmf-chip route add chip.svg --id default -w 0 -w 10 | dmf-chip route list -
found element: `dmf:ChipDesign`
found element: `dmf:TestRoutes`
Add new element: `dmf:TestRoute[@version="0.1.0"][@id="default"]`
Added 2 waypoints.
Wrote test route to: `<stdout>`
Detected Inkscape 0.92+; using 96 pixels per inch
{'version': '0.1.0', 'id': 'default', 'waypoints': ['0', '10']}
```

# Command `route remove`

```sh
Usage: dmf-chip route remove [OPTIONS] CHIP_FILE [OUTPUT]

  Remove routes

Options:
  --id TEXT  id of route to remove. Multiple ids may be specified, e.g.,
             `--id foo --id bar ...`  [required]
  --help     Show this message and exit.
```

### Example

Remove route with `id="default"` from `chip.svg`:

```sh
$ dmf-chip route remove chip.svg --id default | dmf-chip route list -
Removed test route: `default`
Wrote modified chip file to: `<stdout>`
Detected Inkscape 0.92+; using 96 pixels per inch
```

-------------------------------------------------------------------------------

Install
-------

The latest [`dmf-chip` release][1] is available as a
[Conda][2] package from the [`sci-bots`][2] channel.

To install `dmf-chip` in an **activated Conda environment**, run:

    conda install -c sci-bots -c conda-forge dmf-chip

-------------------------------------------------------------------------------

License
-------

This project is licensed under the terms of the [BSD license](/LICENSE.md)

-------------------------------------------------------------------------------

Contributors
------------

 - Christian Fobel ([@sci-bots](https://github.com/sci-bots))


[1]: https://github.com/sci-bots/dmf-chip
[2]: https://anaconda.org/sci-bots/dmf-chip
