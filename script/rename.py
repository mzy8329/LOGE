import os

def rename(dir):
    files = sorted(os.listdir(dir), key=lambda x: x.split("_")[1].split(".")[0])
    for i, file in enumerate(files):
        os.rename(dir+"/"+file, dir+"/"+str(i)+".jpg")
rename("data/imgs")