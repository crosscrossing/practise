#include <iostream>
#include <vector>
#include <deque>
#include <algorithm>
#include <stack>
#include <unordered_map>
#include <queue>
#include <map>
#include <set>

using namespace std;

template<typename T>
void printVector(const vector<T> &i) {
    for (auto &j : i) {
        cout << j << " ";
    }
    cout << std::endl;
}

struct TreeNode {
    int val;
    TreeNode *left;
    TreeNode *right;

    TreeNode() : val(0), left(nullptr), right(nullptr) {}

    TreeNode(int x) : val(x), left(nullptr), right(nullptr) {}

    TreeNode(int x, TreeNode *left, TreeNode *right) : val(x), left(left), right(right) {}

    ~TreeNode() {
        if (left != nullptr) {
            delete left;
        }
        if (right != nullptr) {
            delete right;
        }
    }
};

struct ListNode {
    int val;
    ListNode *next;

    ListNode() : val(0), next(nullptr) {}

    ListNode(int x) : val(x), next(nullptr) {}

    ListNode(int x, ListNode *next) : val(x), next(next) {}
};

class Node {
public:
    int val;
    Node *left;
    Node *right;
    Node *next;

    Node() : val(0), left(NULL), right(NULL), next(NULL) {}

    Node(int _val) : val(_val), left(NULL), right(NULL), next(NULL) {}

    Node(int _val, Node *_left, Node *_right, Node *_next)
            : val(_val), left(_left), right(_right), next(_next) {}
};

struct Status {
    int a, b, c;
};

class Solution {
    //前序遍历
public:
    vector<int> preorderTraversal(TreeNode *root) {
        vector<int> result;
        if (!root) {
            return result;
        }
        stack<TreeNode *> nodeStack;
        nodeStack.push(root);

        while (!nodeStack.empty()) {
            auto node = nodeStack.top();
            nodeStack.pop();
            result.push_back(node->val);

            if (node->right) {
                nodeStack.push(node->right);
            }
            if (node->left) {
                nodeStack.push(node->left);
            }
        }
        return result;
    }

    //中序遍历
    //通过栈遍历
    vector<int> inorderTraversal(TreeNode *root) {
        vector<int> result;
        if (root == nullptr) {
            return result;
        }
        TreeNode *curNode = root; //优先处理
        stack<TreeNode *> nodeStack; //随后处理
        while (curNode != nullptr || !nodeStack.empty()) {
            //优先处理left
            if (curNode != nullptr) {
                nodeStack.push(curNode);
                curNode = curNode->left;
            }
                //随后处理自己
            else if (!nodeStack.empty()) {
                TreeNode *node = nodeStack.top();
                result.push_back(node->val);
                nodeStack.pop();
                // 自己处理完成 下一个处理right
                curNode = node->right;
            }
        }
        return result;
    }

    //中序遍历 // morris遍历
    vector<int> inorderTraversal2(TreeNode *root) {
        vector<int> result;
        if (root == nullptr) {
            return result;
        }
        TreeNode *cur = root, *pre = nullptr;
        while (cur) {
            if (!cur->left) {
                result.push_back(cur->val);
                cur = cur->right;
                continue;
            }
            pre = cur->left;
            while (pre->right && pre->right != cur) {
                pre = pre->right;
            }
            if (!pre->right) {
                pre->right = cur;
                cur = cur->left;
            } else {
                pre->right = nullptr;
                result.push_back(cur->val);
                cur = cur->right;
            }
        }
        return result;
    }

    // 后序
    vector<int> postorderTraversal(TreeNode *root) {
        stack<TreeNode *> postOrderStack;
        if (root) {
            postOrderStack.push(root);
        }
        vector<int> result;
        while (!postOrderStack.empty()) {
            auto node = postOrderStack.top();
            if (node->left == nullptr && node->right == nullptr) {
                postOrderStack.pop();
                result.push_back(node->val);
                delete node;
            } else {
                if (node->right) {
                    postOrderStack.push(node->right);
                    node->right = nullptr;
                }
                if (node->left) {
                    postOrderStack.push(node->left);
                    node->left = nullptr;
                }
            }
        }
        return result;
    }

    // 中后序遍历恢复树
private:
    unordered_map<int, int> inorderMap;

    TreeNode *buildChildren(vector<int> &inorder, vector<int> &postorder, int front, int back, int &post) {
        //构造本节点
        auto *node = new TreeNode(postorder[post]);

        //找自己
        int cur = postorder[post];
        int curInorder = inorderMap[cur];

        while (--post >= 0) {
            //后续从后往前依次找 找孩子
            int child = postorder[post];
            int childInorder = inorderMap[child];
            //范围内没找到 超界了 退回post
            if (childInorder < front) {
                post++;
                break;
            }
            //中序中 孩子在左侧
            if (childInorder < curInorder) {
                node->left = buildChildren(inorder, postorder, front, curInorder - 1, post);
            }
            //中序中 孩子在右侧
            if (childInorder > curInorder) {
                node->right = buildChildren(inorder, postorder, curInorder + 1, back, post);
            }
        }
        return node;
    }

public:
    TreeNode *buildTree(vector<int> &inorder, vector<int> &postorder) {
        if (postorder.empty()) {
            return nullptr;
        }
        for (auto i = 0; i < inorder.size(); i++) {
            inorderMap.insert({inorder[i], i});
        }
        int post = postorder.size() - 1;
        return buildChildren(inorder, postorder, 0, inorder.size() - 1, post);
    }

//相机覆盖二叉树
private:
    Status dfs(TreeNode *root) {
        if (!root) {
            return {INT_MAX / 2, 0, 0};
        }
        auto left = dfs(root->left);
        auto right = dfs(root->right);
        //a withCam ：当前子树 root 有相机，情况下的minCam
        //b noCamWatchByDad：当前子树 root 没有相机，被父亲监控，情况下的minCam
        //c noCamWatchBySon ：当前子树 root 没有相机，被儿子监控，情况下的minCam

        //当前节点 root 放了相机（当前子树的相机数，保底为1）3种情况最小值
        //左右儿子都没有放相机，都被父亲监控
        //minCam(root.left, false, true) + minCam(root.right, false, true)
        //左儿子放了相机，被监控，右儿子没有相机，被父亲监控
        //minCam(root.left, true, true) + minCam(root.right, false, true)
        //左儿子没有相机，被父亲监控，右儿子放了相机，被监控
        //minCam(root.left, false, true) + minCam(root.right, true, true)

        int a = 1 + min(left.b + right.b, min(left.a + right.b, left.b + right.a));

        //当前节点 root 没有相机，但被父亲监控了 4种情况最小值
        //两个儿子都放了相机，被监控
        //左儿子放了相机，被监控，右儿子没有相机，没有被父亲和自己监控，但被自己儿子监控
        //右儿子放了相机，被监控，左儿子没有相机，没有被父亲和自己监控，但被自己儿子监控
        //两个儿子都没有相机，没有被父亲和自己监控，但被自己儿子监控
        int b = min(left.a + right.a, min(left.a + right.c, min(left.c + right.a, left.c + right.c)));

        //当前节点 root 没有相机，也没有被父亲监控，是被儿子监控 3种情况最小值
        //两个儿子都放了相机，被监控
        //左儿子有相机，被监控，右儿子没有相机，没被父亲和自己监控，被自己儿子监控
        //左儿子没有相机，没被父亲和自己监控，被自己儿子监控。右儿子有相机，被监控
        int c = min(left.a + right.a, min(left.a + right.c, left.c + right.a));
        return {a, b, c};
    }

public:
    int minCameraCover(TreeNode *root) {
        auto result = dfs(root);
        return min(result.a, result.c);
    }

//合并二×树
private:
    bool mergeDfs(TreeNode *t1, TreeNode *t2, TreeNode *t) {
        if (!t1 && !t2) {
            delete t;
            return false;
        }
        t->val = (t1 ? t1->val : 0) + (t2 ? t2->val : 0);
        t->left = new TreeNode(0);
        if (!mergeDfs(t1 ? t1->left : nullptr, t2 ? t2->left : nullptr, t->left)) {
            t->left = nullptr;
        }
        t->right = new TreeNode(0);
        if (!mergeDfs(t1 ? t1->right : nullptr, t2 ? t2->right : nullptr, t->right)) {
            t->right = nullptr;
        }
        return true;
    }

public:
    TreeNode *mergeTrees(TreeNode *t1, TreeNode *t2) {
        if (!t1 && !t2) {
            return nullptr;
        }
        TreeNode *t = new TreeNode(0);
        mergeDfs(t1, t2, t);
        return t;
    }

//找排序树里出现频率最多的
private:
    int base, count, maxCount;
    vector<int> answer;

    void update(int x) {
        if (x == base) {
            ++count;
        } else {
            count = 1;
            base = x;
        }
        if (count == maxCount) {
            answer.push_back(base);
        }
        if (count > maxCount) {
            maxCount = count;
            answer = vector<int>{base};
        }
    }

public:
    vector<int> findMode(TreeNode *root) {
        TreeNode *cur = root, *pre = nullptr;
        while (cur) {
            if (!cur->left) {
                update(cur->val);
                cur = cur->right;
                continue;
            }
            pre = cur->left;
            while (pre->right && pre->right != cur) {
                pre = pre->right;
            }
            if (!pre->right) {
                pre->right = cur;
                cur = cur->left;
            } else {
                pre->right = nullptr;
                update(cur->val);
                cur = cur->right;
            }
        }
        return answer;
    }

    // +1
public:
    vector<int> plusOne(vector<int> &digits) {
        int add = 1;
        for (auto it = digits.rbegin(); it < digits.rend() && add; it++) {
            *it += add--;
            if (*it > 9) {
                *it = 0;
                add++;
            }
        }
        if (add > 0) {
            vector<int> newDigits = {1};
            newDigits.insert(newDigits.end(), digits.begin(), digits.end());
            digits = newDigits;
        }
        return digits;
    }

public:
    TreeNode *lowestCommonAncestor(TreeNode *root, TreeNode *p, TreeNode *q) {
//        if (root->val > p->val && root->val < q->val) {
//            return root;
//        }
//        else if (root->val == p->val || root->val == q->val) {
//            return root;
//        }
//        else
        if (root->val > p->val && root->val > q->val) {
            return lowestCommonAncestor(root->left, p, q);
        } else if (root->val < p->val && root->val < q->val) {
            return lowestCommonAncestor(root->right, p, q);
        }
        return root;
    }

private:
    void getPath(TreeNode *node, int sum, vector<int> &path, vector<vector<int>> &result) {
        if (node == nullptr) {
            return;
        }
        sum -= node->val;
        path.push_back(node->val);
        if (sum == 0 && node->left == nullptr && node->right == nullptr) {
            result.push_back(path);
        } else {
            this->getPath(node->left, sum, path, result);
            this->getPath(node->right, sum, path, result);
        }
        path.pop_back();
    }

public:
    vector<vector<int>> pathSum(TreeNode *root, int sum) {
        vector<int> path;
        vector<vector<int>> result;
        this->getPath(root, sum, path, result);
        return result;
    }

// node指向下一个同行节点
private:
    queue<Node *> q;

    void pushChild(Node *next) {
        if (next->left) {
            q.push(next->left);
        }
        if (next->right) {
            q.push(next->right);
        }
    }

public:
    Node *connect(Node *root) {
        if (root) { q.push(root); }
        while (!q.empty()) {
            int size = q.size();
            Node *cur = nullptr;
            Node *last = nullptr;
            for (int i = 0; i < size; i++) {
                cur = q.front();
                q.pop();
                pushChild(cur);
                if (i > 0) {
                    last->next = cur;
                }
                last = cur;
            }
        }
        return root;
    }

public:
    //排序树最小差值
    int getMinimumDifference(TreeNode *root) {
        auto nodes = this->inorderTraversal(root);
        int min = 0x0fffffff;
        for (int i = 0; i < nodes.size() - 1; i++) {
            min = nodes[i + 1] - nodes[i] < min ? nodes[i + 1] - nodes[i] : min;
        }
        return min;
    }


public:
    //链表交换
    ListNode *swapPairs(ListNode *head) {
        auto *begin = new ListNode(0, head);
        auto cur = begin;
        while (cur->next && cur->next->next) {
            auto temp = cur->next;
            cur->next = cur->next->next;
            temp->next = cur->next->next;
            cur->next->next = temp;
            cur = temp;
        }
        return begin->next;
    }

public:
    // 统计每个string都出现的所有char
    vector<string> commonChars(vector<string> &A) {
        map<string, int> countResult;
        for (int i = 0; i < A.size(); ++i) {
            map<string, int> count;
            for (auto &c : A[i]) {
                count[string{c}]++;
            }
            if (i == 0) {
                countResult = count;
            } else {
                for (auto &p : countResult) {
                    countResult[p.first] =
                            count[p.first] < countResult[p.first] ? count[p.first] : countResult[p.first];
                }
            }
        }
        vector<string> result;
        for (auto &p : countResult) {
            for (int i = 0; i < p.second; i++) {
                result.push_back(p.first);
            }
        }
        return result;
    }

public:
    // 求二次方后排序
    vector<int> sortedSquares(vector<int> &A) {
        vector<int> result;
        result.reserve(A.size());
        for (int i = 0, j = A.size() - 1; i <= j;) {
            if (abs(A[i]) < abs(A[j])) {
                result.push_back(i * i);
                i++;
            } else {
                result.push_back(j * j);
                j--;
            }
        }
        sort(result.begin(), result.end());
        return result;

        int size = A.size();
        std::vector<int> ans(A.size());

        int i = 0, j = size - 1;
        while (i <= j) {
            if (abs(A[i]) < abs(A[j])) {
                ans[j - i] = A[j] * A[j];
                j--;
            } else {
                ans[j - i] = A[i] * A[i];
                i++;
            }
        }

        return ans;
    }

    //长按键入 输入重复字母满足条件
public:
    bool isLongPressedName(string name, string typed) {
        auto nameChar = name.rbegin();
        auto typedChar = typed.rbegin();
        char cur = 0, last = 0;
        if (nameChar != name.rend()) {
            cur = *nameChar;
        }
        while (typedChar != typed.rend()) {
            if (*typedChar == cur) {
                ++nameChar;
                last = cur;
                cur = 0;
                if (nameChar != name.rend()) {
                    cur = *nameChar;
                }
                ++typedChar;
            } else if (*typedChar == last) {
                ++typedChar;
            } else {
                return false;
            }
        }
        if (nameChar != name.rend()) {
            return false;
        }
        return true;
    }

public:
    vector<int> partitionLabels(string S) {

    }

    // 找出比当前数更小的有几个
private:
    static bool smaller(pair<int, int> &i, pair<int, int> &j) {
        if (i.second < j.second) {
            return true;
        }
        return false;
    }

public:
    vector<int> smallerNumbersThanCurrent(vector<int> &nums) {
        vector<pair<int, int>> sortedNums;
        for (int i = 0; i < nums.size(); i++) {
            sortedNums.emplace_back(i, nums[i]);

        }
        sort(sortedNums.begin(), sortedNums.end(), smaller);

        vector<int> result(nums.size());
        for (int i = 0; i < nums.size(); i++) {
            int j = i;
            while (j > 0 && sortedNums[j].second == sortedNums[j - 1].second) {
                j--;
            }
            result[sortedNums[i].first] = j;
        }
        return result;
    }

    // 每个数出现的次数不同
public:
    bool uniqueOccurrences(vector<int> &arr) {
        map<int, int> countMap;
        for (auto &i : arr) {
            countMap[i]++;
        }
        vector<int> count;
        for (auto &i : countMap) {
            count.push_back(i.second);
        }
        sort(count.begin(), count.end());
        for (int i = 0; i < count.size() - 1; i++) {
            if (count[i] == count[i + 1]) {
                return false;
            }
        }
        return true;
    }

    // 深度遍历得到的数求和
private:
    int getSum(TreeNode *root, int prev) {
        prev += root->val;
        int sum = 0;
        if (!root->left && !root->right) {
            sum += prev;
        }
        if (root->left) {
            sum += getSum(root->left, prev * 10);
        }
        if (root->right) {
            sum += getSum(root->right, prev * 10);
        }

        return sum;
    }

public:
    int sumNumbers(TreeNode *root) {
        if (!root) {
            return 0;
        }
        return getSum(root, 0);
    }

    // 岛的周长
public:
    int islandPerimeter(vector<vector<int>> &grid) {
        int perimeter = 0;
        int width = grid[0].size(), height = grid.size();
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                if (grid[i][j]) {
                    if (i - 1 < 0 || i - 1 >= height || grid[i - 1][j] == 0) {
                        perimeter++;
                    }
                    if (i + 1 < 0 || i + 1 >= height || grid[i + 1][j] == 0) {
                        perimeter++;
                    }
                    if (j - 1 < 0 || j - 1 >= width || grid[i][j - 1] == 0) {
                        perimeter++;
                    }
                    if (j + 1 < 0 || j + 1 >= width || grid[i][j + 1] == 0) {
                        perimeter++;
                    }
                }
            }
        }
        return perimeter;
    }

public:
    vector<int> intersection(vector<int> &nums1, vector<int> &nums2) {
        set<int> num;
        for (auto &i: nums1) {
            num.insert(i);
        }
        set<int> result;
        for (auto &i: nums2) {
            if (num.find(i) != num.end()) {
                result.insert(i);
            }
        }
        vector<int> res;
        res.reserve(result.size());
        for (auto &i : result) {
            res.push_back(i);
        }
        return res;
    }


public:
    //    /**
    //     * 思路：
    //     *  贪心算法：每个R阵营的参议员禁止下一个离他最近的D阵营的参议员，反之亦然。
    //     * 算法流程：
    //     *  使用两个队列分别保存R阵营和D阵营的参议员索引，
    //     *  在每一轮循环中，比较相邻两个R和D阵营的参议员的索引，
    //     *      保留索引小（min）的，并将该(min + senate.length)添加到该阵营队列尾部
    //     *      去除索引大的，即不添加到末尾。
    //     */
    string predictPartyVictory(string senate) {
        queue<int> rq, dq;
        for (int i = 0; i < senate.size(); i++) {
            switch (senate[i]) {
                case 'R':
                    rq.push(i);
                    break;
                case 'D':
                    dq.push(i);
                    break;
            }
        }

        while (!rq.empty() && !dq.empty()) {
            int rf = rq.front(), df = dq.front();
            rq.pop();
            dq.pop();
            if (rf < df) {
                rq.push(rf + senate.size());
            } else {
                dq.push(df + senate.size());
            }
        }

        if (rq.empty()) {
            return "Dire";
        }
        if (dq.empty()) {
            return "Radiant";
        }
        return "";
    }

// 找出t多的字母
public:
    char findTheDifference(string s, string t) {
        int sSum[123] = {}, tSum[123] = {};
        for (char &c : s) {
            sSum[c]++;
        }
        for (char &c : t) {
            tSum[c]++;
        }

        for (int i = 'a'; i <= 'z'; i++) {
            if (tSum[i] > sSum[i]) {
                return i;
            }
        }
        return 0;
    }

// 动态规划：找上楼最小消耗，一次最多走两步
public:
    int minCostClimbingStairs(vector<int>& cost) {
        vector<int> costSum(cost.size(), 0);
        // init
        costSum[0] = cost[0];
        costSum[1] = cost[1];
        for (int i = 2; i < cost.size(); i++) {
            costSum[i] = cost[i] + (costSum[i-2] > costSum[i-1] ? costSum[i-1] : costSum[i-2]);
        }
        return costSum[cost.size()-2] > costSum[cost.size()-1] ? costSum[cost.size()-1] : costSum[cost.size()-2];
    }

    int climbStairs(int n) {
        if (n == 1) {
            return 1;
        }
        vector<int> ways(n, 0);
        //init
        ways[0] = 1;
        ways[1] = 2;
        for (int i = 2; i < ways.size(); i++) {
            ways[i] = ways[i-1] + ways[i-2];
        }
        return ways[n-1];
    }

//[103] 二叉树的锯齿形层序遍历
public:
    vector<vector<int>> zigzagLevelOrder(TreeNode* root) {
        deque<TreeNode*> inverse, reverse;
        vector<vector<int>> result;

        //init
        inverse.push_back(root);
        for(int i = 0; !inverse.empty() || !reverse.empty(); i++) {
            vector<int> level;
            if (i % 2 == 0) {
                while (!inverse.empty()) {
                    TreeNode* node = inverse.back();
                    level.push_back(node->val);
                    inverse.pop_back();
                    if (node->left) {
                        reverse.push_front(node->left);
                    }
                    if (node->right) {
                        reverse.push_front(node->right);
                    }
                }
            }
            else {
                while (!reverse.empty()) {
                    TreeNode* node = reverse.front();
                    level.push_back(node->val);
                    reverse.pop_front();
                    if (node->right) {
                        inverse.push_back(node->right);
                    }
                    if (node->left) {
                        inverse.push_back(node->left);
                    }
                }
            }
            if (!level.empty()) {
                result.push_back(level);
            }
        }
        return result;
    }
};

int main() {
    Solution s;

    vector<int> aa = {9, 9, 9, 9, 9, 9, 9, 9};
    printVector(s.plusOne(aa));

    auto *t1 = new TreeNode(1);
    t1->left = new TreeNode(3);
    t1->left->left = new TreeNode(5);
    t1->right = new TreeNode(2);
    auto *t2 = new TreeNode(2);
    t2->left = new TreeNode(1);
    t2->left->right = new TreeNode(4);
    t2->right = new TreeNode(3);
    t2->right->right = new TreeNode(7);

    auto t = s.mergeTrees(t1, t2);
    auto tp = s.preorderTraversal(t);
    auto ti = s.inorderTraversal2(t);


    vector<int> inorder{3, 6, 2, 9, 8, 7, 5, 10, 4}, postorder{6, 3, 8, 9, 7, 4, 10, 5, 2};
    auto tree = s.buildTree(inorder, postorder);

    auto *t3 = new TreeNode(6);
    t3->left = new TreeNode(2);
    t3->left->left = new TreeNode(0);
    t3->left->right = new TreeNode(4);
    t3->left->right->left = new TreeNode(3);
    t3->left->right->right = new TreeNode(5);

    auto *p = new TreeNode(0);
    auto *q = new TreeNode(3);
    auto node = s.lowestCommonAncestor(t3, p, q);

    //[5,4,8,11,null,13,4,7,2,null,null,5,1]
    //22

    auto *t4 = new TreeNode(1);
    t4->left = new TreeNode(-2);
    t4->right = new TreeNode(-3);
    t4->left->left = new TreeNode(1);
    t4->left->right = new TreeNode(3);
    t4->right->left = new TreeNode(-2);
    t4->left->left->left = new TreeNode(-1);

    //[5,4,8,11,null,13,4,7,2,null,null,5,1]
    //22

    auto *t5 = new TreeNode(5);
    t5->left = new TreeNode(4);
    t5->right = new TreeNode(8);
    t5->left->left = new TreeNode(11);
    t5->right->left = new TreeNode(13);
    t5->right->right = new TreeNode(4);
    t5->left->left->left = new TreeNode(7);
    t5->left->left->right = new TreeNode(2);
    t5->right->right->left = new TreeNode(5);
    t5->right->right->right = new TreeNode(1);
    auto pathSum = s.pathSum(t5, 22);

    Node *node1 = new Node(1);
    node1->left = new Node(2);
    node1->right = new Node(3);
    node1->left->left = new Node(5);
    node1->right->left = new Node(4);
    s.connect(node1);

    auto ppostorder = s.postorderTraversal(t3);

    auto *t6 = new TreeNode(4);
    t6->left = new TreeNode(2);
    t6->right = new TreeNode(6);
    t6->left->left = new TreeNode(1);
    t6->left->right = new TreeNode(3);
    t6->right->left = new TreeNode(5);
    t6->left->left->left = new TreeNode(0);
    cout << s.getMinimumDifference(t6) << endl;

    auto *linklist = new ListNode(1);
    linklist->next = new ListNode(2);
    linklist->next->next = new ListNode(3);
//    linklist->next->next->next = new  ListNode(4);
    auto newlist = s.swapPairs(linklist);

    vector<string> strings1{"abcc", "cabc", "1"};
    printVector(s.commonChars(strings1));

    cout << s.isLongPressedName("ppyplrza", "pyypllrza") << endl;

    vector<int> nums{6, 6, 6, 7, 7, 7};
    printVector(s.smallerNumbersThanCurrent(nums));

    cout << s.uniqueOccurrences(nums) << endl;

    cout << s.sumNumbers(t6) << endl;

    vector<vector<int>> island{{1, 1}};
    cout << s.islandPerimeter(island) << endl;

    cout << s.predictPartyVictory("RDDRRDDR") << endl;

    cout << s.findTheDifference("abcdabcd", "abcdedbac") << endl;

    vector<int> stairs{1,100,1,2,100,100,7,2,100};
    cout << s.minCostClimbingStairs(stairs) << endl;

    s.zigzagLevelOrder(t6);
    return 0;
}
