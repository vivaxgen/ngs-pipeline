import itertools


class DNACodec:

    # versions
    # 0: 91-char coded, composed of Base-84 Block Scheme and 7-char deletion scheme

    def __init__(self, version: str = "0"):
        """Initializes the codec with an optional default encoding version."""
        self.encoding_version = version

        # 1. Global shared resources (the base 91 printable ASCII alphabet)
        self.ascii_91 = [
            chr(i) for i in range(33, 127) if chr(i) not in ('"', "'", "\\")
        ]

        self.nucs = ["A", "C", "G", "T"]

        # 2. Initialize version-specific data structures
        self._init_v0()
        # Future version initializers can be called conditionally or automatically:
        # if hasattr(self, f"_init_v{self.encoding_version}"):
        #     getattr(self, f"_init_v{self.encoding_version}")()

        # Validate that the requested encoding version exists
        if not hasattr(self, f"_encode_v{self.encoding_version}"):
            raise ValueError(f"Unsupported encoding version: {self.encoding_version}")

    def _init_v0(self):
        """Sets up tables and mappings specific to Version 0 (Base-84 Block Scheme)."""
        triplets = ["".join(p) for p in itertools.product(self.nucs, repeat=3)]
        self.triplet_to_char_v0 = {
            triplet: self.ascii_91[i] for i, triplet in enumerate(triplets)
        }
        self.char_to_triplet_v0 = {v: k for k, v in self.triplet_to_char_v0.items()}

        remainder_options = [
            "".join(p) for p in itertools.product(self.nucs, repeat=1)
        ] + ["".join(p) for p in itertools.product(self.nucs, repeat=2)]
        self.rem_to_char_v0 = {
            rem: self.ascii_91[64 + i] for i, rem in enumerate(remainder_options)
        }
        self.char_to_rem_v0 = {v: k for k, v in self.rem_to_char_v0.items()}

        self.max_del_len_v0 = 91 - 84  # 7 characters available
        self.del_to_char_v0 = {
            i: self.ascii_91[84 + (i - 1)] for i in range(1, self.max_del_len_v0 + 1)
        }
        self.char_to_del_v0 = {v: k for k, v in self.del_to_char_v0.items()}

    def encode(self, dna_string: str, version: str | None = None) -> str:
        """Encodes a DNA string using the instance's default version or an explicitly passed version."""
        version = version or self.encoding_version
        encoder_method = getattr(self, f"_encode_v{version}", None)

        if not encoder_method:
            raise ValueError(f"Unsupported encoding version: {version}")

        return version + encoder_method(dna_string)

    def decode(self, encoded_string: str) -> str:
        """Inspects the version prefix and dynamically dispatches to the correct decoder."""
        if not encoded_string:
            raise ValueError("Empty encoded string.")

        version = encoded_string[0]
        payload = encoded_string[1:]

        decoder_method = getattr(self, f"_decode_v{version}", None)

        if not decoder_method:
            raise ValueError(f"Unknown version prefix: {version}")

        return decoder_method(payload)

    # ==========================================
    # VERSION 0 LOGIC
    # ==========================================
    def _encode_v0(self, dna_string: str) -> str:
        dna_string = dna_string.upper()
        encoded = []
        i = 0
        n = len(dna_string)

        while i < n:
            if dna_string[i] == "-":
                start = i
                while (
                    i < n and dna_string[i] == "-" and (i - start) < self.max_del_len_v0
                ):
                    i += 1
                run_length = i - start
                encoded.append(self.del_to_char_v0[run_length])
            else:
                chunk = []
                while i < n and dna_string[i] != "-" and len(chunk) < 3:
                    chunk.append(dna_string[i])
                    i += 1
                chunk_str = "".join(chunk)

                if len(chunk_str) == 3:
                    encoded.append(self.triplet_to_char_v0[chunk_str])
                else:
                    encoded.append(self.rem_to_char_v0[chunk_str])

        return "".join(encoded)

    def _decode_v0(self, payload: str) -> str:
        decoded = []
        for char in payload:
            if char in self.char_to_triplet_v0:
                decoded.append(self.char_to_triplet_v0[char])
            elif char in self.char_to_rem_v0:
                decoded.append(self.char_to_rem_v0[char])
            elif char in self.char_to_del_v0:
                decoded.append("-" * self.char_to_del_v0[char])
            else:
                raise ValueError(f"Invalid character in v0 payload: {char}")
        return "".join(decoded)


def test_codec():
    codec = DNACodec(version="0")
    test_cases = [
        "ACTGAT------CGAT-AA",
        "A-C-G-T-",
        "ACGTACGTACGT",
        "----",
        "A--C--G--T--",
        "",
    ]

    for original in test_cases:
        encoded = codec.encode(original)
        decoded = codec.decode(encoded)
        assert original == decoded, f"Failed for {original}: got {decoded}"
        print(f"Original: {original}, Encoded: {encoded}, Decoded: {decoded}")


# EOF
