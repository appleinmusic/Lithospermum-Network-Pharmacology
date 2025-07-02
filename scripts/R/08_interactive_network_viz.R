# 09_interactive_network_viz.R

# --- 1. 设置环境和加载包 ---

# 检查并安装所需的包
packages_to_install <- c("visNetwork", "htmlwidgets", "igraph", "dplyr", "scales")
for (pkg in packages_to_install) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "https://cran.r-project.org/")
  }
}

# 加载包
library(visNetwork)
library(igraph)
library(dplyr)
library(htmlwidgets)
library(scales)

cat("--- 环境设置完成，已加载所需R包 ---\n")

# --- 2. 定义文件路径 ---

# 输入文件
ppi_network_path <- "results/network/ppi_network.rds"
node_stats_path <- "results/tables/key_nodes_analysis.csv" # 包含中心性指标的文件
mapping_path <- "results/tables/target_string_mapping.csv" # 包含ID映射的文件

# 输出文件
output_dir <- "results/figures"
output_html_path <- file.path(output_dir, "interactive_ppi_network.html")

# 确保输出目录存在
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

cat(paste("输入PPI网络文件:", ppi_network_path, "\n"))
cat(paste("输入节点统计文件:", node_stats_path, "\n"))
cat(paste("输入ID映射文件:", mapping_path, "\n"))
cat(paste("输出HTML文件:", output_html_path, "\n"))


# --- 3. 加载并准备数据 ---

# 加载igraph对象
if (!file.exists(ppi_network_path)) {
  stop("错误: PPI网络文件不存在! 请先运行 '03_network_construction.R' 脚本。")
}
g <- readRDS(ppi_network_path)

# 加载节点统计数据
if (!file.exists(node_stats_path)) {
  stop("错误: 节点统计文件不存在! 请先运行 '04_network_analysis_and_viz.R' 脚本。")
}
node_stats <- read.csv(node_stats_path)

# 加载ID映射文件
if (!file.exists(mapping_path)) {
  stop("错误: ID映射文件不存在! 请先运行 '03_network_construction.R' 脚本。")
}
id_mapping <- read.csv(mapping_path)

cat("--- 数据加载成功 ---\n")


# --- 4. 准备visNetwork所需的数据格式 ---

# 核心修复：通过映射文件连接指标和网络节点
# 1. 将指标数据 (node_stats) 与 ID映射表 连接
#    (key_nodes_analysis.csv中的'node'列是基因符号, 对应映射文件中的'target_name')
metrics_mapped <- left_join(node_stats, id_mapping, by = c("node" = "target_name"))

# 2. 从igraph对象准备节点列表，其ID是STRING ID
nodes <- data.frame(id = V(g)$name, stringsAsFactors = FALSE)

# 3. 将igraph节点列表与映射后的指标数据连接，这次使用STRING ID进行匹配
nodes <- left_join(nodes, metrics_mapped, by = c("id" = "string_id"))


# --- 诊断 ---
cat("--- 节点数据合并后摘要 ---\n")
print(summary(nodes))

# 检查 'betweenness' 列中是否存在 NA
if(any(is.na(nodes$betweenness))) {
  cat("\n警告: 'betweenness' 列中存在NA值。\n")
}

# 创建颜色映射函数 (更稳健的方式)
pal <- scales::col_numeric("viridis", domain = nodes$betweenness, na.color = "#808080")

# 设置节点的视觉属性
nodes <- nodes %>%
  mutate(
    label = node, # 节点上显示的标签 (基因符号)
    value = degree, # 节点大小，基于度
    title = paste0("<p><b>", node, " (", id, ")</b></p>", # 鼠标悬停时显示的HTML内容
                   "Degree: ", round(degree, 2), "<br>",
                   "Betweenness: ", round(betweenness, 6), "<br>",
                   "Closeness: ", round(closeness, 4)),
    color.background = pal(betweenness), # 颜色，基于betweenness
    color.border = "black",
    color.highlight.background = "orange",
    color.highlight.border = "darkred"
  )

# 准备边数据 (edges data frame)
edges <- igraph::as_data_frame(g, what = "edges")
edges <- edges %>%
    mutate(
        width = combined_score / max(combined_score) * 5 + 1, # 标准化combined_score作为边的宽度
        title = paste("Interaction Score:", combined_score) # 鼠标悬停时显示
    )

cat("--- visNetwork数据格式准备完成 ---\n")
# print(head(nodes))
# print(head(edges))

# --- 5. 创建并配置visNetwork交互式网络图 ---

set.seed(123) # 为了布局的可复现性

vis_net <- visNetwork(nodes, edges, main = "Interactive PPI Network of Lithospermum Targets") %>%
  # 布局设置
  visIgraphLayout(layout = "layout_with_fr") %>%
  # 节点设置
  visNodes(
    shape = "dot",
    size = "value",
    font = list(size = 14)
  ) %>%
  # 边设置
  visEdges(
    smooth = list(enabled = TRUE, type = "dynamic"),
    arrows = "to", # 如果是无向图可以去掉
    color = list(color = "#848484", highlight = "red")
  ) %>%
  # 交互选项
  visOptions(
    highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
    nodesIdSelection = TRUE, # 显示一个下拉菜单来选择节点
    selectedBy = "degree"
  ) %>%
  # 物理引擎设置，用于动态布局
  visPhysics(
    solver = "forceAtlas2Based",
    forceAtlas2Based = list(gravitationalConstant = -50, springLength = 100, springConstant = 0.08)
  ) %>%
  # 添加图例
  visLegend(
    addNodes = list(
      list(label = "High Betweenness", shape = "dot", color = list(background=pal(max(nodes$betweenness, na.rm=TRUE)))),
      list(label = "Low Betweenness", shape = "dot", color = list(background=pal(min(nodes$betweenness, na.rm=TRUE))))
    ),
    useGroups = FALSE,
    main = "Legend"
  ) %>%
  # 开启操作界面
  visInteraction(
    navigationButtons = TRUE, # 显示导航按钮
    dragNodes = TRUE,
    dragView = TRUE,
    zoomView = TRUE,
    keyboard = TRUE
  )

cat("--- 交互式网络图对象创建成功 ---\n")

# --- 6. 保存为HTML文件 ---

saveWidget(vis_net, file = output_html_path, selfcontained = TRUE)

cat(paste0("--- 成功! 交互式网络图已保存到: '", output_html_path, "' ---\n"))
cat("您现在可以在浏览器中打开此文件进行交互式浏览。\n")
