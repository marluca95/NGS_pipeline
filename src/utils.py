from Levenshtein import distance as levenshtein_distance
import difflib
import datetime
import os
import torch


class Logger():
    def __init__(self, directory, file_name: str):
        self.file_name = file_name
        self.directory = directory
        self.created_at = datetime.datetime.now()
        #file_name = f'{file_name}_{self.created_at}'
        with open(f'{directory}/{file_name}', 'w') as f:
            f.truncate(0)
            f.close()
        self.log(f'#---------Logger initiated with name "{file_name}" at {self.created_at}---------#')
    def log(self, text: str):
        device = 'cpu'
        if torch.cuda.is_available():
            device = "cuda:0"
        if(device == 'cpu'):
            with open(f'{self.directory}/{self.file_name}', 'a') as f:
                f.write(text)
                f.write('\n')
                f.close()
        else:
            print(text)

def check_sequence_similarity(ref: str, query: str, max_dist: int, tiled: bool) -> bool:
    """
    Check whether two sequences differ by at most a given edit distance,
    and optionally whether all edits occur within a window of size max_dist if the library was designed in a tiled fashion

    Parameters
    ----------
    ref : str
        Reference sequence.
    query : str
        Query sequence.
    max_dist : int
        Maximum allowed edit distance.
    tiled : bool, optional
        If True, require that all edits occur within a window of length max_dist.

    Returns
    -------
    bool
        True if conditions are met, otherwise False.
    """
    import difflib

    # Quick exit for trivial equality
    if ref == query:
        return True

    dist = levenshtein_distance(ref, query)
    if dist > max_dist:
        return False

    # --- If tiled=True, check edit locality ---
    if tiled:
        matcher = difflib.SequenceMatcher(None, ref, query)
        edits = []
        for tag, i1, i2, j1, j2 in matcher.get_opcodes():
            if tag != 'equal':
                edits.append((i1, i2, j1, j2))
        if not edits:
            return True  # no differences

        # Find smallest window covering all edits
        min_ref = min(i1 for i1, _, _, _ in edits)
        max_ref = max(i2 for _, i2, _, _ in edits)
        if (max_ref - min_ref) > max_dist:
            return False

    return True