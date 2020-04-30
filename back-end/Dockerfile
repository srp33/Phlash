FROM tiangolo/meinheld-gunicorn-flask:python3.7

COPY requirements.txt /app/

RUN python -m pip install --upgrade pip \
 && pip install -r requirements.txt

COPY *.py /app/
