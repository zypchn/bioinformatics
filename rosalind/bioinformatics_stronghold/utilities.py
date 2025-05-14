def readFile(filePath):
    with open(filePath, "r") as f:
        return [l.strip() for l in f.readlines()]