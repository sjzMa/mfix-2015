/*
Code taken from the example on Wikipedia:  http://en.wikipedia.org/wiki/Select_(Unix)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>

#include <sys/select.h>
#include <fcntl.h>
#include <unistd.h>
#include <err.h>
#include <errno.h>
#include <signal.h>

extern void handle_command(char*, char*, ssize_t*);

extern void flush_err_msg_gui(char*);

void die(const char *msg)
{
  perror(msg);
  exit(EXIT_FAILURE);
}

//-----------------------------------------------
// Global variables
//-----------------------------------------------

int cmd_sockfd, cmd_maxfd;
int log_sockfd, log_maxfd;
fd_set cmd_master, cmd_readfds;
fd_set log_master, log_readfds;

struct timeval timeout;

//-----------------------------------------------
// for testing
//-----------------------------------------------
/*
int main(int argc, char **argv) {
  init_socket("8888");

  while(1) {
      printf("DO SOMETHING USEFUL\n");
      sleep(1);
      check_socket();
  }
  return 0;
}*/

//-----------------------------------------------
//     Initialize socket listening on port
//-----------------------------------------------

void init_log_socket(char* port) {

  struct addrinfo *res0, *res, hints;
  int error;
  int on = 1;

  (void)memset(&hints, '\0', sizeof(struct addrinfo));

  hints.ai_family = AF_INET;
  hints.ai_socktype = SOCK_STREAM;
  hints.ai_protocol = IPPROTO_TCP;
  hints.ai_flags = AI_PASSIVE;

  if(0 != (error = getaddrinfo(NULL, port, &hints, &res0)))
    errx(EXIT_FAILURE, "%s", gai_strerror(error));

  for(res = res0; res; res = res->ai_next)
    {
      if(-1 == (log_sockfd = socket(res->ai_family, res->ai_socktype, res->ai_protocol)))
       {
         perror("socket()");
         continue;
       }

      if(-1 == (setsockopt(log_sockfd, SOL_SOCKET, SO_REUSEADDR, (char*)&on, sizeof(int))))
       {
         perror("setsockopt()");
         continue;
       }

      if(-1 == (bind(log_sockfd, res->ai_addr, res->ai_addrlen)))
       {
         perror("bind");
         continue;
       }

      break;

    }

  if(-1 == log_sockfd)
    exit(EXIT_FAILURE);

  freeaddrinfo(res0);

  if(-1 == (listen(log_sockfd, 32)))
    die("listen()");

  if(-1 == (fcntl(log_sockfd, F_SETFD, O_NONBLOCK)))
    die("fcntl()");

  FD_ZERO(&log_master);
  FD_ZERO(&log_readfds);

  FD_SET(log_sockfd, &log_master);

  log_maxfd = log_sockfd;

  timeout.tv_sec = 0;
  timeout.tv_usec = 0;

  signal(SIGPIPE, SIG_IGN);

}

void init_cmd_socket(char* port) {

  struct addrinfo *res0, *res, hints;
  int error;
  int on = 1;

  (void)memset(&hints, '\0', sizeof(struct addrinfo));

  hints.ai_family = AF_INET;
  hints.ai_socktype = SOCK_STREAM;
  hints.ai_protocol = IPPROTO_TCP;
  hints.ai_flags = AI_PASSIVE;

  if(0 != (error = getaddrinfo(NULL, port, &hints, &res0)))
    errx(EXIT_FAILURE, "%s", gai_strerror(error));

  for(res = res0; res; res = res->ai_next)
    {
      if(-1 == (cmd_sockfd = socket(res->ai_family, res->ai_socktype, res->ai_protocol)))
       {
         perror("socket()");
         continue;
       }

      if(-1 == (setsockopt(cmd_sockfd, SOL_SOCKET, SO_REUSEADDR, (char*)&on, sizeof(int))))
       {
         perror("setsockopt()");
         continue;
       }

      if(-1 == (bind(cmd_sockfd, res->ai_addr, res->ai_addrlen)))
       {
         perror("bind");
         continue;
       }

      break;

    }

  if(-1 == cmd_sockfd)
    exit(EXIT_FAILURE);

  freeaddrinfo(res0);

  if(-1 == (listen(cmd_sockfd, 32)))
    die("listen()");

  if(-1 == (fcntl(cmd_sockfd, F_SETFD, O_NONBLOCK)))
    die("fcntl()");

  FD_ZERO(&cmd_master);
  FD_ZERO(&cmd_readfds);

  FD_SET(cmd_sockfd, &cmd_master);

  cmd_maxfd = cmd_sockfd;

  timeout.tv_sec = 0;
  timeout.tv_usec = 0;

  signal(SIGPIPE, SIG_IGN);

}

//-----------------------------------------------
//     Check for for connections on socket
//-----------------------------------------------

void check_sockets() {

  int nready;
  int new;

  memcpy(&cmd_readfds, &cmd_master, sizeof(cmd_master));
  memcpy(&log_readfds, &log_master, sizeof(log_master));

  if(-1 == (nready = select(log_maxfd+1, &log_readfds, NULL, NULL, &timeout)))
    die("select()");

  for(int ii=0; ii<=log_maxfd && nready>0; ii++) {
    if(FD_ISSET(ii, &log_readfds)) {
      nready--;

      if(ii == log_sockfd) {
	if(-1 == (new = accept(log_sockfd, NULL, NULL))) {
	  if(EWOULDBLOCK != errno)
	    perror("accept()");

	  break;
	} else {
	  if(-1 == (fcntl(new, F_SETFD, O_NONBLOCK)))
	    die("fcntl()");

	  FD_SET(new, &log_master);

	  if(log_maxfd < new)
	    log_maxfd = new;
	}
      } else {
	char outbuffer[BUFSIZ];
	flush_err_msg_gui(outbuffer);

	ssize_t nbytes;
	if ((nbytes = send(ii, outbuffer, strlen(outbuffer), 0)) == -1) {
	  // Client probably disconnected
	  perror("send");
	  close(ii);
	  FD_CLR(ii, &log_master);
	}
      }
    }
  }

  if(-1 == (nready = select(cmd_maxfd+1, &cmd_readfds, NULL, NULL, &timeout)))
    die("select()");

  for(int ii=0; ii<=cmd_maxfd && nready>0; ii++) {
    if(FD_ISSET(ii, &cmd_readfds)) {
      nready--;

      if(ii == cmd_sockfd) {
	if(-1 == (new = accept(cmd_sockfd, NULL, NULL))) {
	  if(EWOULDBLOCK != errno)
	    perror("accept()");

	  break;
	} else {
	  if(-1 == (fcntl(new, F_SETFD, O_NONBLOCK)))
	    die("fcntl()");

	  FD_SET(new, &cmd_master);

	  if(cmd_maxfd < new)
	    cmd_maxfd = new;
	}
      } else {
	char buffer[BUFSIZ];
	char outbuffer[BUFSIZ];

	ssize_t nbytes;
	if ((nbytes = recv(ii, buffer, strlen(buffer), 0)) == -1) {
	  if(nbytes <= 0) {
	    if(EWOULDBLOCK != errno)
	      die("recv()");
	    
	    break;
	  }
	}

	handle_command(outbuffer, buffer, &nbytes);

	if ((nbytes = send(ii, outbuffer, strlen(outbuffer), 0)) == -1) {
	  // Client probably disconnected
	  perror("send");
	  close(ii);
	  FD_CLR(ii, &cmd_master);
	}
      }
    }
  }

}
