
import java.lang.Math.*;
import java.util.Random;

public class GraphAlgo extends Object {

    public static void allShortestPathLength(int n, int m, int nodei[],
                         int nodej[], boolean directed[], int weight[], 
                         int root, int mindistance[])
    {
      int i,j,k,n2,large,nodeu,nodev,minlen,temp,minj=0,minv=0;
      int location[] = new int[n+1];
      int distance[][] = new int[n+1][n+1];

      // obtain a large number greater than all edge weights
      large = 1;
      for (k=1; k<=m; k++)
        large += weight[k];

      // set up the distance matrix
      for (i=1; i<=n; i++)
        for (j=1; j<=n; j++)
          distance[i][j] = (i == j) ? 0 : large;
      for (k=1; k<=m; k++) {
        i = nodei[k];
        j = nodej[k];
        if (directed[k])
          distance[i][j] = weight[k];
        else
          distance[i][j] = distance[j][i] = weight[k];
      }

      if (root != 1) {
        // interchange rows 1 and root
        for (i=1; i<=n; i++) {
          temp = distance[1][i];
          distance[1][i] = distance[root][i];
          distance[root][i] = temp;
        }
        // interchange columns 1 and root
        for (i=1; i<=n; i++) {
          temp = distance[i][1];
          distance[i][1] = distance[i][root];
          distance[i][root] = temp;
        }
      }
      nodeu = 1;
      n2 = n + 2;
      for (i=1; i<=n; i++) {
        location[i] = i;
        mindistance[i] = distance[nodeu][i];
      }
      for (i=2; i<=n; i++) {
        k = n2 - i;
        minlen = large;
        for (j=2; j<=k; j++) {
          nodev = location[j];
          temp = mindistance[nodeu] + distance[nodeu][nodev];
          if (temp < mindistance[nodev]) mindistance[nodev] = temp;
          if (minlen > mindistance[nodev]) {
            minlen = mindistance[nodev];
            minv = nodev;
            minj = j;
          }
        }
        nodeu = minv;
        location[minj] = location[k];
      }
      if (root != 1) {
        mindistance[1] = mindistance[root];
        mindistance[root] = 0;
        // interchange rows 1 and root
        for (i=1; i<=n; i++) {
          temp = distance[1][i];
          distance[1][i] = distance[root][i];
          distance[root][i] = temp;
        }
        // interchange columns 1 and root
        for (i=1; i<=n; i++) {
          temp = distance[i][1];
          distance[i][1] = distance[i][root];
          distance[i][root] = temp;
        }
      }
    }


    public static void allPairsShortestPaths(int n, int dist[][], int big,
                                    int startnode, int endnode, int path[])
    {
      int i,j,k,d,num,node;
      int next[][] = new int[n+1][n+1];
      int order[] = new int[n+1];

      // compute the shortest path distance matrix
      for (i=1; i<=n; i++)
        for (j=1; j<=n; j++)
          next[i][j] = i;
      for (i=1; i<=n; i++)
        for (j=1; j<=n; j++)
          if (dist[j][i] < big)
            for (k=1; k<=n; k++)
              if (dist[i][k] < big) {
                d = dist[j][i] + dist[i][k];
                if (d < dist[j][k]) {
                   dist[j][k] = d;
                   next[j][k] = next[i][k];
                }
              }
      // find the shortest path from startnode to endnode
      if (startnode == 0) return;
      j = endnode;
      num = 1;
      order[num] = endnode;
      while (true) {
        node = next[startnode][j];
        num++;
        order[num] = node;
        if (node == startnode) break;
        j = node;
      }
      for (i=1; i<=num; i++)
        path[i] = order[num-i+1];
      path[0] = num;
    }

    public static void breadthFirstSearch(int n, int m, int nodei[],
                          int nodej[], int parent[], int sequence[])
    {
      int i,j,k,enqueue,dequeue,queuelength,p,q,u,v;
      int queue[] = new int[n+1];
      int firstedges[] = new int[n+2];
      int endnode[] = new int[m+1];
      boolean mark[] = new boolean[m+1];
      boolean iterate,found;

      // set up the forward star representation of the graph
      for (j=1; j<=m; j++)
        mark[j] = true;
      firstedges[1] = 0;
      k = 0;
      for (i=1; i<=n; i++) {
        for (j=1; j<=m; j++)
          if (mark[j]) {
            if (nodei[j] == i) {
              k++;
              endnode[k] = nodej[j];
              mark[j] = false;
            }
            else {
              if (nodej[j] == i) {
                k++;
                endnode[k] = nodei[j];
                mark[j] = false;
              }
            }
          }
        firstedges[i+1] = k;
      }
      for (i=1; i<=n; i++) {
        sequence[i] = 0;
        parent[i] = 0;
      }
      k = 0;
      p = 1;
      enqueue = 1;
      dequeue = 1;
      queuelength = enqueue;
      queue[enqueue] = p;
      k++;
      sequence[p] = k;
      parent[p] = 0;
      iterate = true;
      // store all descendants
      while (iterate) {
        for (q=1; q<=n; q++) {
          // check if p and q are adjacent
          if (p < q) {
            u = p;
            v = q;
          }
          else {
            u = q;
            v = p;
          }
          found = false;
          for (i=firstedges[u]+1; i<=firstedges[u+1]; i++)
            if (endnode[i] == v) {
              // p and q are adjacent
              found = true;
              break;
            }
          if (found && sequence[q] == 0) {
            enqueue++;
            if (n < enqueue) enqueue = 1;
            queue[enqueue] = q;
            k++;
            parent[q] = p;
            sequence[q] = k;
          }
        }
        // process all nodes of the same height
        if (enqueue >= dequeue) {
          if (dequeue == queuelength) {
            queuelength = enqueue;
          }
          p = queue[dequeue];
          dequeue++;
          if (n < dequeue) dequeue = 1;
          iterate = true;
          // process other components
        }
        else {
          iterate = false;
          for (i=1; i<=n; i++)
            if (sequence[i] == 0) {
              dequeue = 1;
              enqueue = 1;
              queue[enqueue] = i;
              queuelength = 1;
              k++;
              sequence[i] = k;
              parent[i] = 0;
              p = i;
              iterate = true;
              break;
            }
        }
      }
    }

    public static void connectedComponents(int n, int m, int nodei[], int nodej[],
                                           int component[])
    {
      int edges,i,j,numcomp,p,q,r,typea,typeb,typec,tracka,trackb;
      int compkey,key1,key2,key3,nodeu,nodev;
      int numnodes[] = new int[n+1];
      int aux[] = new int[n+1];
      int index[] = new int[3];

      typec=0;
      index[1] = 1;
      index[2] = 2;
      q = 2;
      for (i=1; i<=n; i++) {
        component[i] = -i;
        numnodes[i] = 1;
        aux[i] = 0;   
      }   
      j = 1;
      edges = m;
      do {
        nodeu = nodei[j];
        nodev = nodej[j];
        key1 = component[nodeu];
        if (key1 < 0) key1 = nodeu;
        key2 = component[nodev];
        if (key2 < 0) key2 = nodev;
        if (key1 == key2) {
          if(j >= edges) {
            edges--;
            break;
          }
          nodei[j] = nodei[edges];
          nodej[j] = nodej[edges];
          nodei[edges] = nodeu;
          nodej[edges] = nodev;
          edges--;
        }
        else {
          if (numnodes[key1] >= numnodes[key2]) {
            key3 = key1;
            key1 = key2;
            key2 = key3;
            typec = -component[key2];
          }
          else
            typec = Math.abs(component[key2]);
          aux[typec] = key1;
          component[key2] = component[key1];
          i = key1;
          do {
            component[i] = key2;
            i = aux[i];
          } while (i != 0);
          numnodes[key2] += numnodes[key1];
          numnodes[key1] = 0;
          j++;
          if (j > edges || j > n) break;
        }
      } while (true);
      numcomp = 0;
      for (i=1; i<=n; i++)
        if (numnodes[i] != 0) {
          numcomp++;
          numnodes[numcomp] = numnodes[i];
          aux[i] = numcomp;
        }
      for (i=1; i<=n; i++) {
        key3 = component[i];
        if (key3 < 0) key3 = i;
        component[i] = aux[key3];
      }      
      if (numcomp == 1) {
        component[0] = numcomp;
        return;
      }
      typeb = numnodes[1];
      numnodes[1] = 1;
      for (i=2; i<=numcomp; i++) {
        typea = numnodes[i];
        numnodes[i] = numnodes[i-1] + typeb - 1;
        typeb = typea;
      }
      for (i=1; i<=edges; i++) {
        typec = nodei[i];
        compkey = component[typec];
        aux[i] = numnodes[compkey];
        numnodes[compkey]++;
      }
      for (i=1; i<=q; i++) {
        typea = index[i];
        do {
          if (typea <= i) break;
          typeb = index[typea];
          index[typea] = -typeb;
          typea = typeb;
        } while (true);
        index[i] = -index[i];
      }
      if (aux[1] >= 0)
        for (j=1; j<=edges; j++) {
          tracka = aux[j];
          do {
            if (tracka <= j) break;
            trackb = aux[tracka];
            aux[tracka] = -trackb;
            tracka = trackb;
          } while (true);
          aux[j] = -aux[j];
        }
      for (i=1; i<=q; i++) {
        typea = -index[i];
        if(typea >= 0) {
          r = 0;
          do {
            typea = index[typea];
            r++;
          } while (typea > 0);
          typea = i;
          for (j=1; j<=edges; j++)
            if (aux[j] <= 0) {
              trackb = j;
              p = r;
              do {
                tracka = trackb;
                key1 = (typea == 1) ? nodei[tracka] : nodej[tracka];
                do {
                  typea = Math.abs(index[typea]);
                  key1 = (typea == 1) ? nodei[tracka] : nodej[tracka];
                  tracka = Math.abs(aux[tracka]);
                  key2 = (typea == 1) ? nodei[tracka] : nodej[tracka];
                  if (typea == 1)
                    nodei[tracka] = key1;
                  else
                    nodej[tracka] = key1;
                  key1 = key2;
                  if (tracka == trackb) {
                    p--;
                    if (typea == i) break; 
                  }
                } while (true);
                trackb = Math.abs(aux[trackb]);
              } while (p != 0);
            }
        }
      }
      for (i=1; i<=q; i++)
        index[i] = Math.abs(index[i]);
      if (aux[1] > 0) {
        component[0] = numcomp;
        return;
      }
      for (j=1; j<=edges; j++)
        aux[j] = Math.abs(aux[j]);
      typea=1;
      for (i=1; i<=numcomp; i++) {
        typeb = numnodes[i];
        numnodes[i] = typeb - typea + 1;
        typea = typeb;
      }
      component[0] = numcomp;
    }

    public static void depthFirstSearch(int n, int m, int nodei[],
                        int nodej[], int parent[], int sequence[])
    {
      int i,j,k,counter,stackindex,p,q,u,v;
      int stack[] = new int[n+1];
      int firstedges[] = new int[n+2];
      int endnode[] = new int[m+1];
      boolean mark[] = new boolean[m+1];
      boolean skip,found;

      // set up the forward star representation of the graph
      for (j=1; j<=m; j++)
        mark[j] = true;
      firstedges[1] = 0;
      k = 0;
      for (i=1; i<=n; i++) {
        for (j=1; j<=m; j++)
          if (mark[j]) {
            if (nodei[j] == i) {
              k++;
              endnode[k] = nodej[j];
              mark[j] = false;
            }
            else {
              if (nodej[j] == i) {
                k++;
                endnode[k] = nodei[j];
                mark[j] = false;
              }
            }
          }
        firstedges[i+1] = k;
      }
      for (i=1; i<=n; i++) {
        sequence[i] = 0;
        parent[i] = 0;
        stack[i] = 0;
      } 
      counter = 0;
      p = 1;
      stackindex = 0;
      // process descendents of node p
      while (true) {
        skip = false;
        counter++;
        parent[p] = 0;
        sequence[p] = counter;
        stackindex++;
        stack[stackindex] = p;
        while (true) {
          skip = false;
          q = 0;
          while (true) {
            q++;
            if (q <= n) {
              // check if p and q are adjacent
              if (p < q) {
                u = p;
                v = q;
              }
              else {
                u = q;
                v = p;
              }
              found = false;
              for (k=firstedges[u]+1; k<=firstedges[u+1]; k++)
                if (endnode[k] == v) {
                  // u and v are adjacent
                  found = true;
                  break;
                }
              if (found && sequence[q] == 0) {
                stackindex++;
                stack[stackindex] = q;
                parent[q] = p;
                counter++;
                sequence[q] = counter;
                p = q;
                if (counter == n) return;
                break;
              }
            }
            else {
              // back up
              stackindex--;
              if (stackindex > 0) {
                q = p;
                p = stack[stackindex];
              }
              else {
                skip = true;
                break;
              }
            }
          }
          if (skip) break;
        }
        // process the next component
        stackindex = 0;
        skip = false;
        for (k=1; k<=n; k++)
          if (sequence[k] == 0) {
            p = k;
            skip = true;
            break;
          }
        if (!skip) break;
      }
    }

    public static void heapsort(int n, int x[], boolean ascending)
    {
      /* sort array elements x[1], x[2],..., x[n] in the order of
       *  increasing (ascending=true) or decreasing (ascending=false)
       */

      int elm,h,i,index,k,temp;

      if (n <= 1) return;
      // initially nodes n/2 to 1 are leaves in the heap
      for (i=n/2; i>=1; i--) {
        elm = x[i];
        index = i;
        while (true) {
          k = 2 * index;
          if (n < k) break;
          if (k + 1 <= n) {
            if (ascending) {
              if (x[k] < x[k+1]) k++;
            }
            else {
              if (x[k+1] < x[k]) k++;
            }
          }
          if (ascending) {
            if (x[k] <= elm) break;
          }
          else {
            if (elm <= x[k]) break;
          }
          x[index] = x[k];
          index = k;
        }
        x[index] = elm;
      }
      // swap x[1] with x[n]
      temp = x[1];
      x[1] = x[n];
      x[n] = temp;
      // repeat delete the root from the heap 
      for (h=n-1; h>=2; h--) {
        // restore the heap structure of x[1] through x[h]
        for (i=h/2; i>=1; i--) {
          elm = x[i];
          index = i;
          while (true) {
            k = 2 * index;
            if (h < k) break;
            if (k + 1 <= h) {
              if (ascending) {
                if (x[k] < x[k+1]) k++;
              }
              else {
                if (x[k+1] < x[k]) k++;
              }
            }
            if (ascending) {
              if (x[k] <= elm) break;
            }
            else {
              if (elm <= x[k]) break;
            }
            x[index] = x[k];
            index = k;
          }
          x[index] = elm;
        }
        // swap x[1] and x[h]
        temp = x[1];
        x[1] = x[h];
        x[h] = temp;
      }
    }

    public static void minimumSpanningTreeKruskal(int n, int m, int nodei[], 
                  int nodej[], int weight[], int treearc1[], int treearc2[])
    {
      int i,index,index1,index2,index3,nodeu,nodev,nodew,len,nedge,treearc;
      int halfm,numarc,nedge2=0;
      int predecessor[] = new int[n+1];

      for (i=1; i<=n; i++)
        predecessor[i] = -1;
      // initialize the heap structure
      i = m / 2;
      while (i > 0) {
        index1 = i;
        halfm = m / 2;
        while (index1 <= halfm) {
          index = index1 + index1;
          index2 = ((index < m) && (weight[index + 1] < weight[index])) ?
                   index + 1 : index;
          if (weight[index2] < weight[index1]) {
            nodeu = nodei[index1];
            nodev = nodej[index1];
            len   = weight[index1];
            nodei[index1] = nodei[index2];
            nodej[index1] = nodej[index2];
            weight[index1] = weight[index2];
            nodei[index2] = nodeu;
            nodej[index2] = nodev;
            weight[index2] = len;
            index1 = index2;
          }
          else
            index1 = m;
        }
        i--;
      }
      nedge = m;
      treearc = 0;
      numarc = 0;
      while ((treearc < n-1) && (numarc < m)) {
        // examine the next edge
        numarc++;
        nodeu = nodei[1];
        nodev = nodej[1];
        nodew = nodeu;
        // check if nodeu and nodev are in the same component
        while (predecessor[nodew] > 0) {
          nodew = predecessor[nodew];
        }
        index1 = nodew;
        nodew = nodev;
        while (predecessor[nodew] > 0) {
          nodew = predecessor[nodew];
        }
        index2 = nodew;
        if (index1 != index2) {
          // include nodeu and nodev in the minimum spanning tree
          index3 = predecessor[index1] + predecessor[index2];
          if (predecessor[index1] > predecessor[index2]) {
            predecessor[index1] = index2;
            predecessor[index2] = index3;
          }
          else {
            predecessor[index2] = index1;
            predecessor[index1] = index3;
          }
          treearc++;
          treearc1[treearc] = nodeu;
          treearc2[treearc] = nodev;
        }
        // restore the heap structure
        nodei[1] = nodei[nedge];
        nodej[1] = nodej[nedge];
        weight[1] = weight[nedge];
        nedge--;
        index1 = 1;
        nedge2 = nedge / 2;
        while (index1 <= nedge2) {
          index = index1 + index1;
          index2 = ((index < nedge) && (weight[index + 1] < weight[index])) ?
                   index + 1 : index;
          if (weight[index2] < weight[index1]) {
            nodeu = nodei[index1];
            nodev = nodej[index1];
            len   = weight[index1];
            nodei[index1] = nodei[index2];
            nodej[index1] = nodej[index2];
            weight[index1] = weight[index2];
            nodei[index2] = nodeu;
            nodej[index2] = nodev;
            weight[index2] = len;
            index1 = index2;
          }
          else
            index1 = nedge;
        }
      }
      treearc1[0] = treearc;
    }


    public static void minimumSpanningTreePrim(int n, int dist[][], int tree[])
    {
      int i,j,n1,d,mindist,node,k=0;

      n1 = n - 1;
      for (i=1; i<=n1; i++)
        tree[i] = -n;
      tree[n] = 0;
      for (i=1; i<=n1; i++) {
        mindist = Integer.MAX_VALUE;
        for (j=1; j<=n1; j++) {
          node = tree[j];
          if (node <= 0) {
            d = dist[-node][j];
            if (d < mindist) {
              mindist = d;
              k = j;
            }
          }
        }
        tree[k] = -tree[k];
        for (j=1; j<=n1; j++) {
          node = tree[j];
          if (node <= 0)
            if (dist[j][k] < dist[j][-node]) tree[j] = -k;
        }
      }
    }

    public static void nodeColoring(int n, int m, int nodei[],
                                    int nodej[], int color[])
    {
      int i,j,k,loop,currentnode,newc,ncolor,paint,index,nodek,up,low;
      int degree[] = new int[n+1];
      int choices[] = new int[n+1];
      int maxcolornum[] = new int[n+1];
      int currentcolor[] = new int[n+1];
      int feasiblecolor[] = new int[n+1];
      int firstedges[] = new int[n+2];
      int endnode[] = new int[m+m+1];
      int availc[][] = new int[n+1][n+1];
      boolean more;

      // set up the forward star representation of the graph
      for (i=1; i<=n; i++)
        degree[i] = 0;
      k = 0;
      for (i=1; i<=n; i++) {
        firstedges[i] = k + 1;
        for (j=1; j<=m; j++)
          if (nodei[j] == i) {
            k++;
            endnode[k] = nodej[j];
            degree[i]++;
          }
          else {
            if (nodej[j] == i) {
              k++;
              endnode[k] = nodei[j];
              degree[i]++;
            }
          }
      }
      firstedges[n+1] = k + 1;
      for (i=1; i<=n; i++) {
        feasiblecolor[i] = degree[i] + 1;
        if (feasiblecolor[i] > i) feasiblecolor[i] = i;
        choices[i] = feasiblecolor[i];
        loop = feasiblecolor[i];
        for (j=1; j<=loop; j++)
          availc[i][j] = n;
        k = feasiblecolor[i] + 1;
        if (k <= n)
          for (j=k; j<=n; j++)
            availc[i][j] = 0;
      }
      currentnode = 1;
      // color currentnode
      newc = 1;
      ncolor = n;
      paint = 0;
      more = true;
      do {
        if (more) {
          index = choices[currentnode];
          if (index > paint + 1) index = paint + 1;
          while ((availc[currentnode][newc] < currentnode) && (newc <= index))
            newc++;
          // currentnode has the color 'newc'
          if (newc == index + 1)
            more = false;
          else {
            if (currentnode == n) {
              // a new coloring is found
              currentcolor[currentnode] = newc;
              for (i=1; i<=n; i++)
                color[i] = currentcolor[i];
              if (newc > paint) paint++;
              ncolor = paint;
              if (ncolor > 2) {
              // backtrack to the first node of color 'ncolor'
                index = 1;
                while (color[index] != ncolor)
                  index++;
                j = n;
                while (j >= index) {
                  currentnode--;
                  newc = currentcolor[currentnode];
                  paint = maxcolornum[currentnode];
                  low = firstedges[currentnode];
                  up = firstedges[currentnode + 1];
                  if (up > low) {
                    up--;
                    for (k=low; k<=up; k++) {
                      nodek = endnode[k];
                      if (nodek > currentnode)
                        if (availc[nodek][newc] == currentnode) {
                          availc[nodek][newc] = n;
                          feasiblecolor[nodek]++;
                        }
                    }
                  }
                  newc++;
                  more = false;
                  j--;
                }
                paint = ncolor - 1;
                for (i=1; i<=n; i++) {
                  loop = choices[i];
                  if (loop > paint) {
                    k = paint + 1;
                    for (j=k; j<=loop; j++)
                      if (availc[i][j] == n) feasiblecolor[i]--;
                    choices[i] = paint;
                  }
                }
              }
            }
            else {
              // currentnode is less than n
              low = firstedges[currentnode];
              up = firstedges[currentnode + 1];
              if (up > low) {
                up--;
                k = low;
                while ((k <= up) && more) {
                  nodek = endnode[k];
                  if (nodek > currentnode)
                    more = !((feasiblecolor[nodek] == 1) &&
                           (availc[nodek][newc] >= currentnode));
                  k++;
                }
              }
              if (more) {
                currentcolor[currentnode] = newc;
                maxcolornum[currentnode] = paint;
                if (newc > paint) paint++;
                low = firstedges[currentnode];
                up = firstedges[currentnode + 1];
                if (up > low) {
                  up--;
                  for (k=low; k<=up; k++) {
                    nodek = endnode[k];
                    if (nodek > currentnode)
                      if (availc[nodek][newc] >= currentnode) {
                        availc[nodek][newc] = currentnode;
                        feasiblecolor[nodek]--;
                      }
                  }
                }
                currentnode++;
                newc = 1;
              }
              else
                newc++;
            }
          }
        }
        else {  
          more = true;
          if ((newc > choices[currentnode]) || (newc > paint + 1)) {
            currentnode--;
            newc = currentcolor[currentnode];
            paint = maxcolornum[currentnode];
            low = firstedges[currentnode];
            up = firstedges[currentnode + 1];
            if (up > low) {
              up--;
              for (k=low; k<=up; k++) {
                nodek = endnode[k];
                if (nodek > currentnode)
                  if (availc[nodek][newc] == currentnode) {
                    availc[nodek][newc] = n;
                    feasiblecolor[nodek]++;
                  }
              }
            }
            newc++;
            more = false;
          }
        }
      } while ((currentnode != 1) && (ncolor != 2));
      color[0] = ncolor;
    }

    public static void stronglyConnectedComponents(int n, int m, int nodei[],
                                                int nodej[], int component[])
    {
      int i,j,k,series,stackpointer,numcompoents,p,q,r;
      int backedge[] = new int[n+1];
      int parent[] = new int[n+1];
      int sequence[] = new int[n+1];
      int stack[] = new int[n+1];
      int firstedges[] = new int[n+2];
      int endnode[] = new int[m+1];
      boolean next[] = new boolean[n+1];
      boolean trace[] = new boolean[n+1];
      boolean fresh[] = new boolean[m+1];
      boolean skip,found;

      // set up the forward star representation of the graph
      firstedges[1] = 0;
      k = 0;
      for (i=1; i<=n; i++) {
        for (j=1; j<=m; j++)
          if (nodei[j] == i) {
            k++;
            endnode[k] = nodej[j];
          }
        firstedges[i+1] = k;
      }
      for (j=1; j<=m; j++)
        fresh[j] = true;
      // initialize
      for (i=1; i<=n; i++) {
        component[i] = 0;
        parent[i] = 0;
        sequence[i] = 0;
        backedge[i] = 0;
        next[i] = false;
        trace[i] = false;
      }
      series = 0;
      stackpointer = 0;
      numcompoents = 0;
      // choose an unprocessed node not in the stack
      while (true) {
        p = 0;
        while (true) {
          p++;
          if (n < p) {
            component[0] = numcompoents;
            return;
          }
          if (!trace[p]) break;
        }
        series++;
        sequence[p] = series;
        backedge[p] = series;
        trace[p] = true;
        stackpointer++;
        stack[stackpointer] = p;
        next[p] = true;
        while (true) {
          skip = false;
          for (q=1; q<=n; q++) {
            // find an unprocessed edge (p,q)
            found = false;
            for (i=firstedges[p]+1; i<=firstedges[p+1]; i++)
              if ((endnode[i] == q) && fresh[i]) {
                // mark the edge as processed
                fresh[i] = false;
                found = true;
                break;
              }
            if (found) {
              if (!trace[q]) {
                series++;
                sequence[q] = series;
                backedge[q] = series;
                parent[q] = p;
                trace[q] = true;
                stackpointer++;
                stack[stackpointer] = q;
                next[q] = true;
                p = q;
              }
              else {
                if (trace[q]) {
                  if (sequence[q] < sequence[p] && next[q]) {
                    backedge[p] = (backedge[p] < sequence[q]) ? 
                                   backedge[p] : sequence[q];
                  }
                }
              }
              skip = true;
              break;
            }  
          }
          if (skip) continue;
          if (backedge[p] == sequence[p]) {
            numcompoents++;
            while (true) {
              r = stack[stackpointer];
              stackpointer--;
              next[r] = false;
              component[r] = numcompoents;
              if (r == p) break;
            }
          }
          if (parent[p] != 0) {
            backedge[parent[p]] = (backedge[parent[p]] < backedge[p]) ? 
                                   backedge[parent[p]] : backedge[p];
            p = parent[p];
          }
          else
            break;
        }
      }
    }
}
