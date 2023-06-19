def parse_pyfasta_int(length: int):
    width = len(str(length))
    return [str(number).zfill(width) for number in range(length)]
