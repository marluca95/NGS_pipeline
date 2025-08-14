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