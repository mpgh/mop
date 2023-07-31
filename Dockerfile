FROM python:3.10
LABEL maintainer="etibachelet@gmail.com"

# the exposed port must match the deployment.yaml containerPort value
EXPOSE 80
ENTRYPOINT [ "/usr/local/bin/gunicorn", "mop.wsgi", "-b", "0.0.0.0:80", "--access-logfile", "-", "--error-logfile", "-", "-k", "gevent", "--timeout", "300", "--workers", "2"]

WORKDIR /mop

COPY requirements.txt /mop
RUN pip install --no-cache-dir -r /mop/requirements.txt

COPY . /mop
