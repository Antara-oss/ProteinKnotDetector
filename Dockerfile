FROM python:3.11

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    HOME=/home/user

RUN useradd -m -u 1000 user && \
    mkdir -p /home/user/app && \
    chown -R user:user /home/user

WORKDIR /home/user/app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY --chown=user:user . .

USER user

EXPOSE 7860

CMD ["gunicorn", "-b", "0.0.0.0:7860", "knot_web.wsgi:application", "--timeout", "120", "--workers", "1"]
