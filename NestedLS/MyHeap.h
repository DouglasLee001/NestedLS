#ifndef _MY_HEAP_H
#define _MY_HEAP_H
#include <iostream>
#include <string>
#include <functional>

class MyHeap
{
public:
    MyHeap(long capacity) : cmpCallBack(nullptr)
    {
        _capacity = capacity;
        pos_in_array = new long[_capacity];
        is_in_heap = new long[_capacity];
        array = new long[_capacity];
        std::fill_n(is_in_heap, _capacity, 0);
        count = 0;
    }
    ~MyHeap()
    {
        delete[] pos_in_array;
        pos_in_array = nullptr;
        delete[] array;
        array = nullptr;
        delete[] is_in_heap;
        is_in_heap = nullptr;
    }
    void clear()
    {
        count = 0;
        std::fill_n(is_in_heap, _capacity, 0);
    }
    void registerCmpCallback(std::function<bool(int, int)> fun)
    {
        cmpCallBack = fun;
    }
    void registerScoreCallback(std::function<int(int)> fun)
    {
        scoreCallBack = fun;
    }
    long getPos(long v)
    {
        return pos_in_array[v];
    }
    long getVAtPos(long pos)
    {
        return array[pos];
    }
    long size() const
    {
        return count;
    }
    long capacity() const
    {
        return _capacity;
    }
    void swap(long a, long b)
    {
        long t = array[a];
        array[a] = array[b];
        pos_in_array[array[a]] = a;

        array[b] = t;
        pos_in_array[array[b]] = b;
    }
    bool is_leaf(long pos) const
    {
        return (pos >= count / 2) && (pos < count);
    }
    bool isInHeap(long v)
    {
        return is_in_heap[v] == 1 && pos_in_array[v] < count;
    }
    static long left_child(long pos)
    {
        return 2 * pos + 1;
    }
    static long right_child(long pos)
    {
        return 2 * pos + 2;
    }
    static long parent(long pos)
    {
        return (pos - 1) / 2;
    }
    void shiftDown(long pos)
    {
        while (!is_leaf(pos))
        {
            long l = left_child(pos);
            long r = right_child(pos);
            if (l >= count)
            {
                return;
            }
            if (r < count && cmpCallBack(array[r], array[l]))
            {
                l = r;
            }
            if (cmpCallBack(array[pos], array[l]))
            {
                return;
            }
            swap(pos, l);
            pos = l;
        }
    }
    void insert(long v)
    {
        if (is_in_heap[v] == 1)
        {
            return;
        }
        long cur = count++;
        array[cur] = v;
        pos_in_array[v] = cur;
        is_in_heap[v] = 1;
        while (cur != 0 && cmpCallBack(array[cur], array[parent(cur)]))
        {
            swap(cur, parent(cur));
            cur = parent(cur);
        }
    }
    void modify(long v)
    {
        if (is_in_heap[v] == 0)
        {
            return;
        }
        long pos = pos_in_array[v];
        if (pos >= count)
        {
            return;
        }
        remove(pos);
        insert(v);
    }
    long remove_first()
    {
        if (count != 0)
        {
            long v = array[0];
            is_in_heap[v] = 0;
            swap(0, --count);
            shiftDown(0);
            return v;
        }
        return -1;
    }
    long remove(long pos)
    {
        if (pos == count - 1)
        {
            count--;
        }
        else
        {
            swap(pos, --count);
            while (pos != 0 && cmpCallBack(array[pos], array[parent(pos)]))
            {
                swap(pos, parent(pos));
                pos = parent(pos);
            }
            if (count != 0)
            {
                shiftDown(pos);
            }
        }
        long remove_v = array[count];
        is_in_heap[remove_v] = 0;
        return remove_v;
    }

private:
    std::function<bool(int, int)> cmpCallBack;
    std::function<int(int)> scoreCallBack;
    long *pos_in_array;
    long *is_in_heap;
    long *array;
    long _capacity;
    long count;
};

#endif
