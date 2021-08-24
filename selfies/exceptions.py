class SMILESParserError(ValueError):

    def __init__(self, smiles, reason="N/A", idx=-1):
        self.smiles = smiles
        self.idx = idx
        self.reason = reason

    def __str__(self):
        err_msg = "\n" \
                  "\tSMILES: {smiles}\n" \
                  "\t        {pointer}\n" \
                  "\tIndex:  {index}\n" \
                  "\tReason: {reason}"

        return err_msg.format(
            smiles=self.smiles,
            pointer=(" " * self.idx + "^"),
            index=self.idx,
            reason=self.reason
        )


class EncoderError(Exception):
    pass


class DecoderError(Exception):
    pass
