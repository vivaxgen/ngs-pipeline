import signal


# Define a handler function that raises a TimeoutError
def timeout_handler(signum, frame):
    raise TimeoutError("Input timed out!")


# Set the alarm signal handler
signal.signal(signal.SIGALRM, timeout_handler)


def timed_input(prompt, timeout, default=None):
    try:
        # Start the timer
        signal.alarm(timeout)
        # Get input from user
        user_input = input(prompt)
        # Cancel the alarm if input is received
        signal.alarm(0)
        return user_input
    except TimeoutError:
        print("\n[Timed out]")
        return default


def test():
    # Usage example with a 5-second timeout
    if (result := timed_input("Enter something in 5 seconds: ", 5)) is None:
        pass
    else:
        print(result)


# EOF
