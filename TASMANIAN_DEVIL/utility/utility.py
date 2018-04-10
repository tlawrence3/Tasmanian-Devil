#exportPickle function from EXAMO package
def exportPickle(obj, fileName, mode = 'wb', protocol = -1):
    import pickle
    f = open(fileName, mode)
    pickle.dump(obj, f, protocol = -1)
    f.close()
