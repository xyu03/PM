

/*-------------------------------------------------------
 | RFWStack: stack of elements of type T. Upper bound
 |           n on maximum size must be known in advance.
 |           Initialization done in O(n) time, all other
 |           operations in constant time.
 |
 | author: Microsoft Corporation
 *------------------------------------------------------*/

#ifndef RFW_STACK_H
#define RFW_STACK_H
#include <cassert>

template <class T, bool debug=false> class RFWStack {
	private:
		T *stack;
		int top, size;

	public:
		inline bool isFull() {return (top==size);}

		inline void push (T i) {
			if (debug) assert (top<size);
			stack[++top] = i;
		}

		inline T peekTop() {
			return stack[top];
		}

		inline T pop() {
			if (debug) assert (!isEmpty());
			return (stack[top--]);
		};

		inline bool isEmpty() {return (top==0);};

		inline int getNElements() {return top;}

		// peek element at position i; first allowed position is 1
		inline T peek (int p) {
			if (debug) assert (p>0 && p<=top);
			return (stack[p]);
		}

		RFWStack (int s) {
			if (debug) assert (s>0);
			top = 0;
   			size = s;
			stack = new T [size+1];
		};

		inline void reset() {top = 0;}

		~RFWStack () {delete [] stack;};
};

#endif
