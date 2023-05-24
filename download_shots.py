import urllib.request
import zipfile
import os

def download_shot(shotid):
    target_dir = os.path.join('shots', shotid)
    if not os.path.exists(target_dir):
        download_url = f'http://service.iris.edu/ph5ws/dataselect/1/query?reqtype=SHOT&format=segy1&net=8J&sta=*&cha=DPZ&starttime=2019-05-01T00:00:00&endtime=2019-05-31T00:00:00&shotid={shotid}&length=60&nodata=404'
        path, _ = urllib.request.urlretrieve(download_url, reporthook=print)
        with zipfile.ZipFile(path, 'r') as zip_ref:
            print(f'extracting files from {path} into {target_dir}')
            zip_ref.extractall(target_dir)
        urllib.request.urlcleanup()
