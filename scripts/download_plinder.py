import argparse

from plinder.core.utils.gcs import download_dataset


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Download plinder dataset ...')

    parser.add_argument('-b', '--bucket_name', type=str,
                        default="plinder",
                        help='Bucket name if different from \"plinder\"')
    parser.add_argument('-r', '--release', type=str,
                        default="2024-06",
                        help='Release date')
    parser.add_argument('-i', '--iteration', type=str,
                        default="v2",
                        help='Iteration')
    parser.add_argument('-d', '--specific_dirs', nargs='+',
                        default=[],
                        help='List of specific directory to download')
    parser.add_argument('-u', '--unpack', action='store_true',
                        help='unpack zips')
    parser.add_argument('-s', '--skip_download', action='store_true',
                        help='Skip download if the goal is to just unpack already downloaded data')


    args = parser.parse_args()
    print(f"bucket_name: {args.bucket_name}",
            f"release: {args.release}",
            f"iteration: {args.iteration}",
            f"specific_dirs: {args.specific_dirs}",
            f"unpack: {args.unpack}"
            f"skip_download: {args.skip_download}")
    download_dataset(
            bucket_name= args.bucket_name,
            release = args.release,
            iteration = args.iteration,
            specific_dirs = args.specific_dirs,
            unpack = args.unpack,
            skip_download = args.skip_download)
