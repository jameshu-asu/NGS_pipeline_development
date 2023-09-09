#!/usr/bin/env python
import schedule
import time
import subprocess

def job_that_executes_once():
    subprocess.call(['python', 'HAdV_respiratory_Monitor_v3.py'])
    return schedule.CancelJob

schedule.every().day.at('23:00').do(job_that_executes_once)

while True:
    schedule.run_pending()
    time.sleep(1)






