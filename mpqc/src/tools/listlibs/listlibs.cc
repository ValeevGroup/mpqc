
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/wait.h>

////////////////////////////////////////////////////////////////////////

class node {
  public:
    char * val_;
    node * n;
    node * p;

    node(const char *s);
    ~node() { if (val_) free(val_); }

    void print();
  };

node::node(const char *s) : val_(0), n(0), p(0)
{
  // remove quotes and newlines from s if they exist
  if (s[0] == '"') {
    val_ = strdup(&s[1]);
    val_[strlen(val_)-2] = '\0';
    }
  else {
    val_ = strdup(s);
    val_[strlen(val_)-1] = '\0';
    }
  }

void node::print()
{
  if (val_) printf("%s ",val_);
  if (n) n->print();
}

//////////////////////////////////////////////////////////////////////

class list {
  private:
    node * head_;
    node * cur_;
  
  public:
    list() : head_(0), cur_(0) {}

    void add(node*);
    void del_cur();

    node * current() { return cur_; }
    int cur_is_on_list_after_cur();

    void rewind() { cur_ = head_; }
    void ff() { while(cur_->n) cur_ = cur_->n; }
    void backspace() { cur_ = cur_->p; }
    void space() { cur_ = cur_->n; }
    void print();

  };

void list::print()
{
  if (head_) head_->print();
}

void list::add(node*nd)
{
  if (!head_) {
    head_ = nd;
    cur_ = head_;
    }
  else {
    ff();
    cur_->n = nd;
    nd->p = cur_;
    cur_ = cur_->n;
    }
  }

void list::del_cur()
{
  node *t;
  if (cur_->n && cur_->p) { // somewhere in the middle
    cur_->p->n = cur_->n;
    cur_->n->p = cur_->p;
    t = cur_;
    cur_ = cur_->n;
    delete t;
    }
  else if (!cur_->n && !cur_->p) { // all done
    delete cur_;
    head_=cur_=0;
    }
  else if (!cur_->n) { // the end
    cur_->p->n = 0;
    t = cur_;
    cur_ = cur_->p;
    delete t;
    }
  else if(!cur_->p) { // the head
    cur_->n->p = 0;
    t = cur_;
    cur_ = cur_->n;
    head_ = cur_;
    delete t;
    }
  else {
    fprintf(stderr,"what? where are you?\n");
    exit (1);
    }
  }

int list::cur_is_on_list_after_cur()
{
  node *l;

  // start from cur->n, and test to see if the current val already exists

  for (l=cur_->n; l; l=l->n)
    if (!strcmp(cur_->val_,l->val_)) return 1;

  return 0;
  }

////////////////////////////////////////////////////////////////////////

int
run_child(int argc, char *argv[])
{
  int fd[2];

  if (pipe(fd) < 0) {
    fprintf(stderr,"pipe() failed\n");
    exit(1);
    }

  int pid;
  if ((pid=fork()) < 0) {
    fprintf(stderr,"fork() failed\n");
    exit(1);
    }

  if (!pid) { /* child */

    // we don't need the read end of the FIFO
    close(fd[0]);

    // want stdout to go into the FIFO
    dup2(fd[1],1);
    
    // create the arguments for execvp
    char **args = new char*[argc+2];
    args[0] = strdup("gcc");
    args[1] = strdup("-E");

    for (int i=0; i < argc-1 ; i++)
      args[i+2] = strdup(argv[i+1]);

    args[argc+1] = 0;

    // and do it
    execvp("gcc",args);

    // this no good, you fix
    exit(1);
    }

  close(fd[1]);
  return(fd[0]);

  }

void
read_junk(int fd)
{
  char line[1024];
  FILE *str = fdopen(fd,"r");

  list l;

  while(!feof(str)) {
    fgets(line,1023,str);
    if (line[0] != '#' && line[0] != '\n')
      l.add(new node(line));
    }

  //l.print();
  //printf("\n\n");

  l.ff();
  l.backspace();

  while (l.current()) {
    if (l.cur_is_on_list_after_cur()) l.del_cur();
    l.backspace();
    }

  l.print();
  printf("\n\n");
  }
    
main(int argc, char *argv[])
{
  int fd = run_child(argc,argv);

  read_junk(fd);

  /* clean up zombies */

  while(wait(&argc) > -1) ;

  exit(0);
  }
