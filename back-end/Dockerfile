FROM tiangolo/meinheld-gunicorn-flask:python3.7

ADD requirements.txt /app/

RUN python -m pip install --upgrade pip \
 && pip install -r requirements.txt

ADD *.py /app/
ADD genemark_suite_linux_64.tar.gz /
ADD gm_key_64.gz /

RUN gunzip /gm_key_64.gz \
 && mv /gm_key_64 /.gm_key

USER 1000
