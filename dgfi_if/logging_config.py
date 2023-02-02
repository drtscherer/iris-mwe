import logging, sys, os
from datetime import datetime

log_date = datetime.strftime(datetime.now(),'%Y%m%d')


# if 'Logs' not in os.listdir('.'):
#     os.mkdir('Logs')
# if 'dgfi_if_logs' not in os.listdir('Logs/'):
#     os.mkdir('Logs/dgfi_if_logs')


logging.basicConfig(
    format="%(asctime)s %(name)s %(levelname)s : %(message)s",
    datefmt='%Y-%m-%d %H:%M:%S',
    handlers=[
        # logging.FileHandler(f"Logs/dgfi_if_logs/{log_date}.log"),
        logging.StreamHandler(sys.stdout)
    ]
)

level = logging.INFO