#!/usr/bin/env python3
"""Nicolas Gampierakis (2019). Downloads external weather data.
"""

import socket
import sys
import time
import urllib.error
import urllib.request


def progress_bar(packet_number, packet_size, total_size):
    """Custom reporthook for displaying download progress bar. Called as a
    default keyword argument by download_database().

    Args:
        packet_number (int): number of data packets transferred from server
            to client
        packet_size (int): size of data packet
        total_size (int): total database file size reported by Content-Length
            header
    """

    # Create timer to calculate download speed
    global start_tick
    if packet_number == 0:
        start_tick = time.time()
        return
    elapsed_time = time.time() - start_tick
    # Track total downloaded data
    packet_progress = packet_number * packet_size
    # Use min() to avoid percentage overshoot
    percentage_downloaded = min(packet_number * packet_size * 100 /
                                total_size, 100)
    bandwidth = int(packet_progress / (1024 * elapsed_time))
    # Using write instead of print to avoid automatic \n
    sys.stdout.write("\rDownloading weather data ...%d%%, %d MB, %d KB/s" %
                     (percentage_downloaded, packet_progress / (1024 * 1024),
                      bandwidth))
    # Flushes buffer to terminal to avoid lag in output
    sys.stdout.flush()


def download_database(server_path, client_path):
    """Downloads and extracts database into client directory.

    Args:
        client_path (pathlib.PurePath): the file path to the client's
        downloaded database
        server_path (str): the url path to the server database
    """

    # Timeout value recommended by NN/g before user loses attention
    if sys.version_info.major < 3 or sys.version_info.minor < 5:
        raise Exception("Python 3.5 or greater is required to run the update "
                        "package")
    socket.setdefaulttimeout(10)
    data = urllib.request.urlretrieve(server_path, client_path,
                                      reporthook=progress_bar)
    print("\nWeather database successfully downloaded.")
    return data
