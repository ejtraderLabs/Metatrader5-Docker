FROM alpine:3.12 AS st-builder

RUN apk add --no-cache make gcc git freetype-dev \
            fontconfig-dev musl-dev xproto libx11-dev \
            libxft-dev libxext-dev
RUN git clone https://github.com/DenisKramer/st.git /work
WORKDIR /work
RUN make

FROM alpine:3.12 AS xdummy-builder

RUN apk add --no-cache make gcc freetype-dev \
            fontconfig-dev musl-dev xproto libx11-dev \
            libxft-dev libxext-dev
RUN apk add --no-cache libressl3.1-libcrypto libressl3.1-libssl linux-headers --no-cache --repository http://dl-3.alpinelinux.org/alpine/edge/main/
RUN apk add x11vnc --no-cache --repository http://dl-3.alpinelinux.org/alpine/edge/community/
RUN Xdummy -install

# ----------------------------------------------------------------------------

FROM ejtrader/pyzmq:dev

USER root
ENV WINEPREFIX=/root/.wine
ENV WINEARCH=win64
ENV DISPLAY :0
ENV USER=root
ENV PASSWORD=root


# Basic init and admin tools
RUN apk --no-cache add supervisor sudo wget \
    && echo "$USER:$PASSWORD" | /usr/sbin/chpasswd \
    && rm -rf /apk /tmp/* /var/cache/apk/*

# Install X11 server and dummy device
RUN apk add --no-cache xorg-server xf86-video-dummy \
    && apk add libressl3.1-libcrypto --no-cache --repository http://dl-3.alpinelinux.org/alpine/edge/main/ \
    && apk add libressl3.1-libssl --no-cache --repository http://dl-3.alpinelinux.org/alpine/edge/main/ \
    && apk add x11vnc --no-cache --repository http://dl-3.alpinelinux.org/alpine/edge/community/ \
    && rm -rf /apk /tmp/* /var/cache/apk/*
COPY --from=xdummy-builder /usr/bin/Xdummy.so /usr/bin/Xdummy.so
COPY assets/xorg.conf /etc/X11/xorg.conf
COPY assets/xorg.conf.d /etc/X11/xorg.conf.d

# Configure init
COPY assets/supervisord.conf /etc/supervisord.conf

# Openbox window manager
RUN apk --no-cache add openbox  \
    && rm -rf /apk /tmp/* /var/cache/apk/*
COPY assets/openbox/mayday/mayday-arc /usr/share/themes/mayday-arc
COPY assets/openbox/mayday/mayday-arc-dark /usr/share/themes/mayday-arc-dark
COPY assets/openbox/mayday/mayday-grey /usr/share/themes/mayday-grey
COPY assets/openbox/mayday/mayday-plane /usr/share/themes/mayday-plane
COPY assets/openbox/mayday/thesis /usr/share/themes/thesis
COPY assets/openbox/rc.xml /etc/xdg/openbox/rc.xml
COPY assets/openbox/menu.xml /etc/xdg/openbox/menu.xml
COPY Metatrader /root/Metatrader
# Login Manager
RUN apk --no-cache add slim consolekit \
    && rm -rf /apk /tmp/* /var/cache/apk/*
RUN /usr/bin/dbus-uuidgen --ensure=/etc/machine-id
COPY assets/slim/slim.conf /etc/slim.conf
COPY assets/slim/alpinelinux /usr/share/slim/themes/alpinelinux

# A decent system font
RUN apk add --no-cache font-noto \
    && rm -rf /apk /tmp/* /var/cache/apk/*
COPY assets/fonts.conf /etc/fonts/fonts.conf



# st  as terminal
RUN apk add --no-cache freetype fontconfig xproto libx11 libxft libxext ncurses \
    && rm -rf /apk /tmp/* /var/cache/apk/*
COPY --from=st-builder /work/st /usr/bin/st
COPY --from=st-builder /work/st.info /etc/st/st.info
RUN tic -sx /etc/st/st.info

# Some other resources
RUN apk add --no-cache xset \
    && rm -rf /apk /tmp/* /var/cache/apk/*
COPY assets/xinit/Xresources /etc/X11/Xresources
COPY assets/xinit/xinitrc.d /etc/X11/xinit/xinitrc.d

COPY assets/x11vnc-session.sh /root/x11vnc-session.sh
COPY assets/start.sh /root/start.sh


RUN apk update && apk add samba-winbind wine && ln -s /usr/bin/wine64 /usr/bin/wine

# Download Winetricks
RUN wget https://raw.githubusercontent.com/Winetricks/winetricks/master/src/winetricks && chmod +x winetricks && mv winetricks /usr/bin/winetricks
#RUN winetricks -q vcrun2015 dotnet40
# Download Mono
#RUN wget -P /mono http://dl.winehq.org/wine/wine-mono/5.1.1/wine-mono-5.1.1-x86.msi
# Install Mono Runtime for .NET Applications
#RUN wine msiexec /i /mono/wine-mono-5.1.1-x86.msi
# Download gecko
#RUN wget -P /gecko http://dl.winehq.org/wine/wine-gecko/2.47.1/wine-gecko-2.47.1-x86_64.msi
#RUN wine msiexec /i /gecko/wine-gecko-2.47.1-x86_64.msi
#RUN wget -P /vrun https://download.microsoft.com/download/0/6/4/064F84EA-D1DB-4EAA-9A5C-CC2F0FF6A638/vc_redist.x64.exe



WORKDIR /$HOME/
EXPOSE 5900 15555 15556 15557 15558
CMD ["/usr/bin/supervisord","-c","/etc/supervisord.conf"]


