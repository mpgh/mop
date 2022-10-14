import logging
import sys

def start_log():

    # Initialize logged with DEBUG-level output
    log = logging.getLogger('mop_log')
    log.setLevel(logging.DEBUG)

    # Format the default output string:
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    # Add a stream handler to output to STDOUT because this will be accessible
    # through the Kubernetes logging output
    ch = logging.StreamHandler(stream=sys.stdout)
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(formatter)
    log.addHandler(ch)

    return log

def stop_log(log):
    for handler in log.handlers:
        handler.close()
        log.removeFilter(handler)
