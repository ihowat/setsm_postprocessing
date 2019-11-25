# Erik Husby; Polar Geospatial Center, University of Minnesota; 2018


import pyproj


def main():
    parser = argparse.ArgumentParser(
        description="Run the pyproj.transform command and print the results.")
    parser.add_argument('projstr_in', type=str,
        help="Projection of input points.")
    parser.add_argument('projstr_out', type=str,
        help="Projection of output points.")
    parser.add_argument('coordinates', nargs='*', type=float,
        help="All n coordinates for first dimension followed by all n coordinates for second dimension.")

    args = parser.parse_args()

    num_coords = len(args.coordinates)
    if num_coords % 2 != 0:
        parser.error("Odd number of input coordinates were supplied")

    midpoint = int(num_coords / 2)

    x_in = args.coordinates[:midpoint]
    y_in = args.coordinates[midpoint:]

    x_out, y_out = pyproj.transform(pyproj.Proj(args.projstr_in), pyproj.Proj(args.projstr_out), x_in, y_in)

    print(x_out + y_out)

    sys.exit(0)



if __name__ == '__main__':
    main()
